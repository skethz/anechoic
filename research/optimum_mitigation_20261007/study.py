"""Driver for PROTOCOL.md (optimum_mitigation_20261007).
Usage: python3 study.py {a1|a2|c|b|hold|hold_b|b48}
Outputs: results/<stage>.json (+ .npz with per-trial arrays), logs via stdout."""
import hashlib
import itertools
import json
import math
import sys
import time
from pathlib import Path

import numpy as np

import sca_fast as F

ROOT = Path(__file__).resolve().parent
RES = ROOT / 'results'
RES.mkdir(exist_ok=True)
BASE = 20261007
E = 12
THR = (33000, 33200, 33300, 33320, 33337)
FREQS = {'V80_250MHz': 250e6, 'ASIC_1.0GHz': 1.0e9, 'ASIC_1.25GHz': 1.25e9}
C_RS = 886.0
J, SUMW, J8 = F.load()
J2 = (2 * J8).astype(np.int8)


def log(*a):
    print(time.strftime('%H:%M:%S'), *a, flush=True)


def file_hashes():
    return {p: hashlib.sha256((ROOT / p).read_bytes()).hexdigest() for p in ('PROTOCOL.md', 'sca_fast.py', 'study.py')}


def key(cfg):
    return json.dumps(cfg, sort_keys=True)


def S_total(cfg):
    return int(cfg['S']) + int(cfg.get('fin_L', 0))


# ------------------------------------------------------------------------------------------------ statistics / TTS
def tts12(tR_s, p):
    """Eq. 9 with P_R = 1 - (1 - p)^12; returns seconds."""
    if p <= 0:
        return math.inf
    return F.tts_ms(tR_s, 1.0 - (1.0 - p) ** E)


def tts1(t_s, p):
    if p <= 0:
        return math.inf
    return F.tts_ms(t_s, p)


def summarize(cfg, cut, flips, Stot, best=None, greedy=None, extra=None):
    cut = np.asarray(cut); n = len(cut)
    cyc = F.cycles(flips, Stot)
    tR = F.round_cycles(cyc, E)
    r = dict(cfg=cfg, n=n, S_total=Stot, mean_cut=float(cut.mean()), sd_cut=float(cut.std()), max_cut=float(cut.max()),
             mean_flips=float(np.mean(flips)), trial_cycles=float(cyc.mean()), round_cycles_E12=tR)
    for thr in THR:
        r[f'k{thr}'] = int((cut >= thr).sum())
        if best is not None:
            r[f'kbest{thr}'] = int((np.asarray(best) >= thr).sum())
        if greedy is not None:
            r[f'kgreedy{thr}'] = int((np.asarray(greedy) >= thr).sum())
    if greedy is not None:
        r['mean_cut_greedy'] = float(np.mean(greedy))
    r['tts'] = tts_block(r['k33337'], n, cyc.mean(), tR)
    if extra:
        r.update(extra)
    return r


def tts_block(k, n, trial_cyc, round_cyc):
    p = k / n; lo, hi = F.wilson(k, n)
    out = dict(p=p, p_lo=lo, p_hi=hi)
    for name, f in FREQS.items():
        tR = round_cyc / f; t1 = trial_cyc / f
        out[name] = dict(t_trial_ms=1e3 * t1, t_round_ms=1e3 * tR,
                         tts12_ms=1e3 * tts12(tR, p), tts12_ms_lo=1e3 * tts12(tR, hi), tts12_ms_hi=1e3 * tts12(tR, lo),
                         tts1_ms=1e3 * tts1(t1, p))
    return out


def score(r, thr, n_key='n'):
    """12-engine V80 TTS99 [ms] at threshold thr with p at its Wilson lower bound."""
    lo, _ = F.wilson(r[f'k{thr}'], r[n_key])
    return 1e3 * tts12(r['round_cycles_E12'] / FREQS['V80_250MHz'], lo)


def run_and_summarize(cfg, B, seed, greedy=True, extra=None):
    t = time.time()
    out = F.run_cfg(J, J2, SUMW, cfg, B, seed, greedy=greedy)
    r = summarize(cfg, out['cut'], out['flips'], S_total(cfg), out['best'], out.get('cut_greedy'), extra)
    r['sec'] = time.time() - t
    return r, out


def dump(name, obj):
    (RES / f'{name}.json').write_text(json.dumps(obj, indent=1) + '\n')


def load(name):
    return json.loads((RES / f'{name}.json').read_text())


# ------------------------------------------------------------------------------------------------ (a) stage 1 grid
def grid_a():
    g = []
    t1s = lambda q: (5.0, (q + 1) / 2, q / 2)
    for S, T0, q, lam in itertools.product((2000, 3600, 8000, 16000), (10.0, 12.0, 15.0), (6.0, 8.0), (0.9, 1.05)):
        for T1 in t1s(q):
            g.append(dict(family='onsager', S=S, T0=T0, q=q, lam=lam, ramp=True, T1=T1))
    for S, T0, q, kap in itertools.product((2000, 3600, 8000, 16000), (12.0, 15.0, 20.0), (6.0, 8.0), (1.5, 2.0)):
        for T1 in t1s(q):
            g.append(dict(family='tecT', S=S, T0=T0, q=q, kappa=kap, ramp=True, T1=T1))
    for S, T0, q in itertools.product((3600, 8000, 16000), (30.0, 40.0), (4.0, 6.0, 8.0)):
        for T1 in t1s(q):
            g.append(dict(family='plain', S=S, T0=T0, q=q, T1=T1))
    return g


REFERENCE = dict(family='onsager', S=8000, T0=12.0, q=8.0, lam=1.05, ramp=True, T1=5.0)


def stage_a1():
    g = grid_a(); log('stage a1:', len(g), 'configurations x 256 runs', file_hashes())
    rows = []
    for i, cfg in enumerate(g):
        r, _ = run_and_summarize(cfg, 256, [BASE, 1, i])
        r['idx'] = i; rows.append(r)
        log('A1', i, key(cfg), f"mean {r['mean_cut']:.1f} k300 {r['k33300']} k320 {r['k33320']} k337 {r['k33337']} "
            f"kg337 {r['kgreedy33337']} kb337 {r['kbest33337']} flips {r['mean_flips']:.0f} ({r['sec']:.1f}s)")
        if i % 20 == 0:
            dump('a1_screen', dict(hashes=file_hashes(), rows=rows))
    dump('a1_screen', dict(hashes=file_hashes(), rows=rows))


def a1_candidates(rows):
    cands = []
    groups = {}
    for r in rows:
        groups.setdefault((r['cfg']['family'], r['cfg']['S']), []).append(r)
    for (fam, S), rs in sorted(groups.items()):
        rs = sorted(rs, key=lambda r: (score(r, 33300), -r['mean_cut']))
        cands += [r['cfg'] for r in rs[:2]]
    return cands


def stage_a2():
    rows = load('a1_screen')['rows']
    cands = a1_candidates(rows)
    log('stage a2:', len(cands), 'candidates x 2048 runs', file_hashes())
    out = []
    for i, cfg in enumerate(cands):
        r, _ = run_and_summarize(cfg, 2048, [BASE, 2, i])
        r['idx'] = i; out.append(r)
        log('A2', i, key(cfg), f"mean {r['mean_cut']:.1f} k300 {r['k33300']} k320 {r['k33320']} k337 {r['k33337']} "
            f"kg337 {r['kgreedy33337']} score337 {score(r, 33337):.1f} ms ({r['sec']:.1f}s)")
        dump('a2_confirm', dict(hashes=file_hashes(), rows=out))
    sel = select_a(out)
    dump('a2_confirm', dict(hashes=file_hashes(), rows=out, selection=sel))
    log('SELECTION', json.dumps(sel))


def select_a(rows):
    sel = {}
    for fam in ('onsager', 'tecT', 'plain'):
        rs = [r for r in rows if r['cfg']['family'] == fam]
        for cls, ok in (('overall', lambda r: True), ('engine4096', lambda r: r['S_total'] <= 4096)):
            pool = [r for r in rs if ok(r)]
            if not pool:
                continue
            for thr in (33337, 33320, 33300):
                if any(r[f'k{thr}'] > 0 for r in pool):
                    break
            best = min(pool, key=lambda r: (score(r, thr), -r['mean_cut']))
            sel[f'{fam}_{cls}'] = dict(cfg=best['cfg'], ranked_at=thr, score_ms=score(best, thr))
    return sel


def heldout_configs():
    sel = load('a2_confirm')['selection']
    cfgs = []
    for k in ('onsager_overall', 'onsager_engine4096', 'tecT_overall', 'tecT_engine4096', 'plain_overall'):
        if k in sel and sel[k]['cfg'] not in cfgs:
            cfgs.append(sel[k]['cfg'])
    if REFERENCE not in cfgs:
        cfgs.append(REFERENCE)
    return sel, cfgs


# ------------------------------------------------------------------------------------------------ (c) finishing
def fin_options(cfg):
    opts = []
    for L in (50, 100, 200, 400):
        for qf, Tf in ((8.0, 4.5), (6.0, 3.5), (4.0, 2.5), (2.0, 1.5), (1.0, 1.0), (12.0, 6.5), (cfg['q'], cfg['T1'])):
            o = dict(fin_L=L, fin_q=qf, fin_T=Tf)
            if o not in opts:
                opts.append(o)
    return opts


def engine_class(cfg):
    return 'engine4096' if int(cfg['S']) <= 4096 else 'extended'


def valid_options(cfg):
    return [o for o in fin_options(cfg) if engine_class(cfg) == 'extended' or int(cfg['S']) + o['fin_L'] <= 4096]


def anneal_state(cfg, B, seed):
    tabs = F.tables(cfg)
    s0, rs = F.make_batch(J, B, seed)
    st = F.State(J, None, s0, rs)
    F.advance(J2, st, tabs, 0, int(cfg['S']))
    return st


def copy_state(st):
    c = F.State.__new__(F.State)
    for a in ('s', 'sp', 'h', 'nlin', 'flips', 'best2', 'rng'):
        setattr(c, a, getattr(st, a).copy())
    return c


def greedy_cuts(st):
    s = st.s.copy(); h = st.h.copy()
    g = F.greedy_descent(J2, s, h)
    return (SUMW + 0.5 * F.sh_sum(s, h)) / 2.0, g


def finish_eval(cfg, st, options):
    """Apply each finishing option to copies of the same final anneal states. Returns list of (option, arrays)."""
    res = []
    base_cut = F.cuts(st, SUMW); base_g, base_gf = greedy_cuts(st)
    res.append((None, dict(cut=base_cut, flips=st.flips.copy(), best=F.best_cuts(st, SUMW), greedy=base_g, gflips=base_gf)))
    for o in options:
        c2 = dict(cfg, **o)
        tabs = F.tables(c2)
        x = copy_state(st)
        F.advance(J2, x, tabs, int(cfg['S']), int(cfg['S']) + o['fin_L'])
        g, gf = greedy_cuts(x)
        res.append((o, dict(cut=F.cuts(x, SUMW), flips=x.flips.copy(), best=F.best_cuts(x, SUMW), greedy=g, gflips=gf)))
    return res


def stage_c():
    sel, cfgs = heldout_configs()
    log('stage c:', len(cfgs), 'configurations x 4096 pilot anneals', file_hashes())
    out = {}
    for i, cfg in enumerate(cfgs):
        t = time.time()
        rows = []
        for ch in range(2):
            st = anneal_state(cfg, 2048, [BASE, 3, i, ch])
            for j, (o, a) in enumerate(finish_eval(cfg, st, valid_options(cfg))):
                if ch == 0:
                    rows.append(dict(opt=o, cut=[a['cut']], flips=[a['flips']], best=[a['best']], greedy=[a['greedy']],
                                     gflips=[a['gflips']]))
                else:
                    for kk in ('cut', 'flips', 'best', 'greedy', 'gflips'):
                        rows[j][kk].append(a[kk])
        summ = []
        for rw in rows:
            o = rw['opt']; c2 = dict(cfg, **(o or {}))
            cat = {kk: np.concatenate(rw[kk]) for kk in ('cut', 'flips', 'best', 'greedy', 'gflips')}
            r = summarize(c2, cat['cut'], cat['flips'], S_total(c2), cat['best'], cat['greedy'],
                          dict(opt=o, mean_greedy_flips=float(cat['gflips'].mean())))
            summ.append(r)
        opts = [r for r in summ if r['opt'] is not None]
        best = min(opts, key=lambda r: (score(r, 33337), -r['mean_cut']))
        out[key(cfg)] = dict(cfg=cfg, rows=summ, selected=best['opt'], selected_score_ms=score(best, 33337),
                             none_score_ms=score(summ[0], 33337))
        log('C', i, key(cfg), f"none k337 {summ[0]['k33337']} kg {summ[0]['kgreedy33337']} | selected {best['opt']} "
            f"k337 {best['k33337']} score {score(best, 33337):.1f} vs none {score(summ[0], 33337):.1f} ms ({time.time() - t:.0f}s)")
        for r in summ:
            log('   ', r['opt'], f"k300 {r['k33300']} k320 {r['k33320']} k337 {r['k33337']} kg337 {r['kgreedy33337']} "
                f"flips {r['mean_flips']:.0f} score {score(r, 33337):.1f}")
        dump('c_pilot', dict(hashes=file_hashes(), configs=out))


# ------------------------------------------------------------------------------------------------ (a)+(c) held-out
def stage_hold():
    sel, cfgs = heldout_configs()
    cp = load('c_pilot')['configs']
    log('stage hold:', len(cfgs), 'configurations x 8192 held-out runs', file_hashes())
    out = {}
    for i, cfg in enumerate(cfgs):
        t = time.time()
        opts = valid_options(cfg)
        acc = None
        for ch in range(4):
            st = anneal_state(cfg, 2048, [BASE, 9, ch])
            res = finish_eval(cfg, st, opts)
            if acc is None:
                acc = [(o, {kk: [v] for kk, v in a.items()}) for o, a in res]
            else:
                for j, (o, a) in enumerate(res):
                    for kk, v in a.items():
                        acc[j][1][kk].append(v)
        summ = []
        for o, a in acc:
            c2 = dict(cfg, **(o or {}))
            cat = {kk: np.concatenate(v) for kk, v in a.items()}
            r = summarize(c2, cat['cut'], cat['flips'], S_total(c2), cat['best'], cat['greedy'],
                          dict(opt=o, mean_greedy_flips=float(cat['gflips'].mean())))
            summ.append((r, cat))
        sel_opt = cp[key(cfg)]['selected'] if key(cfg) in cp else None
        base_r, base_c = summ[0]
        rec = dict(cfg=cfg, no_finish=base_r, selected_option=sel_opt, options=[r for r, _ in summ[1:]])
        for r, cat in summ[1:]:
            if r['opt'] == sel_opt:
                both = int(((cat['cut'] >= 33337) & (base_c['cut'] >= 33337)).sum())
                gain = int(((cat['cut'] >= 33337) & (base_c['cut'] < 33337)).sum())
                loss = int(((cat['cut'] < 33337) & (base_c['cut'] >= 33337)).sum())
                rec['selected'] = r
                rec['paired'] = dict(both=both, finish_only=gain, anneal_only=loss)
        extra = {}
        for r, c in summ[1:]:
            if sel_opt is not None and r['opt'] == sel_opt:
                extra = dict(cut_sel=c['cut'], flips_sel=c['flips'], greedy_sel=c['greedy'])
        np.savez_compressed(RES / f'hold_{i}.npz', cut=base_c['cut'], flips=base_c['flips'], best=base_c['best'],
                            greedy=base_c['greedy'], **extra)
        out[key(cfg)] = rec
        b = base_r
        msg = (f"no-finish k337 {b['k33337']}/{b['n']} kg {b['kgreedy33337']} kbest {b['kbest33337']} "
               f"TTS12 V80 {b['tts']['V80_250MHz']['tts12_ms']:.1f} ms")
        if 'selected' in rec:
            s_ = rec['selected']
            msg += (f" | finish {sel_opt} k337 {s_['k33337']} TTS12 V80 {s_['tts']['V80_250MHz']['tts12_ms']:.1f} ms "
                    f"paired {rec['paired']}")
        log('HOLD', i, key(cfg), msg, f'({time.time() - t:.0f}s)')
        dump('holdout_ac', dict(hashes=file_hashes(), selection=sel, configs=out))


# ------------------------------------------------------------------------------------------------ (b) population annealing
def pa_configs(base):
    cs = []
    for S in (2000, 4000, 8000):
        cs.append(dict(base, S=S, pa='IND', K=1, alpha=0.0))
        for K, al in itertools.product((4, 8, 16), (0.125, 0.25, 0.5, 1.0, 2.0)):
            cs.append(dict(base, S=S, pa='gibbs', K=K, alpha=al))
        cs.append(dict(base, S=S, pa='pca', K=8, alpha=1.0))
    return cs


def logl(s, h, q, T):
    z = (q + s.astype(np.float64) * h.astype(np.float64)) / T
    return np.logaddexp(0.0, -z).sum(1)


def run_pa(cfg, P, Rr, seed, rs_seed):
    """P populations x Rr replicas. Returns per-population arrays."""
    S = int(cfg['S']); K = int(cfg['K'])
    tabs = F.tables(cfg); T = tabs[0]
    beta = lambda t: 1.0 / T[min(t, S - 1)]
    B = P * Rr
    s0, rs = F.make_batch(J, B, seed)
    st = F.State(J, None, s0, rs)
    g = np.random.default_rng(rs_seed)
    bounds = [int(round(k * S / K)) for k in range(K + 1)]
    stage_cyc = np.zeros(P); n_copied = np.zeros(P); n_unique = np.zeros(P); events = 0
    for k in range(K):
        f0 = st.flips.copy()
        F.advance(J2, st, tabs, bounds[k], bounds[k + 1])
        L = bounds[k + 1] - bounds[k]
        stage_cyc += (0.1424 * (st.flips - f0) + 18.09 * L).reshape(P, Rr).max(1)
        if k < K - 1 and cfg['pa'] != 'IND':
            tb, tn = bounds[k + 1], bounds[k + 2]
            Eng = -0.5 * F.sh_sum(st.s, st.h).astype(np.float64)
            db = beta(tn) - beta(tb)
            lw = -cfg['alpha'] * db * Eng
            if cfg['pa'] == 'pca':
                q = float(tabs[1][tb])
                lw = lw + logl(st.s, st.h, q, T[min(tn, S - 1)]) - logl(st.s, st.h, q, T[tb])
            lw = lw.reshape(P, Rr)
            w = np.exp(lw - lw.max(1, keepdims=True)); w /= w.sum(1, keepdims=True)
            cw = np.cumsum(w, 1); cw[:, -1] = 1.0
            pos = (g.random(P)[:, None] + np.arange(Rr)[None, :]) / Rr
            idx = np.minimum((pos[:, :, None] > cw[:, None, :]).sum(-1), Rr - 1)
            gi = (idx + (np.arange(P) * Rr)[:, None]).ravel()
            n_copied += (gi != np.arange(B)).reshape(P, Rr).sum(1)
            srt = np.sort(idx, 1)
            n_unique += 1 + (np.diff(srt, axis=1) != 0).sum(1)   # distinct parents (surviving lineages) per population
            st.s = st.s[gi]; st.sp = st.sp[gi]; st.h = st.h[gi]; st.nlin = st.nlin[gi]
            events += 1
    cut = F.cuts(st, SUMW).reshape(P, Rr)
    best = F.best_cuts(st, SUMW).reshape(P, Rr)
    gc, _ = greedy_cuts(st)
    tot_cyc = (0.1424 * st.flips + 18.09 * S).reshape(P, Rr).max(1)
    return dict(cut=cut, best=best, greedy=gc.reshape(P, Rr), stage_cyc=stage_cyc, ind_cyc=tot_cyc,
                events=events, n_copied=n_copied, n_unique=n_unique, flips=st.flips.reshape(P, Rr))


def pa_summary(cfg, a, Rr, engines=12):
    P = a['cut'].shape[0]
    mux = max(1, Rr // engines)   # replicas time-multiplexed per engine
    r = dict(cfg=cfg, R=Rr, populations=P, events=a['events'], mean_copied_per_event=float(a['n_copied'].mean() / max(1, a['events'])),
             single_replica_k33337=int((a['cut'] >= 33337).sum()), single_replica_n=int(a['cut'].size),
             mean_final_cut=float(a['cut'].mean()), mean_pop_max=float(a['cut'].max(1).mean()))
    if 'n_unique' in a and a['events']:
        r['mean_distinct_parents_per_event'] = float(a['n_unique'].mean() / a['events'])
    for thr in THR:
        r[f'k{thr}'] = int((a['cut'].max(1) >= thr).sum())
        r[f'kbest{thr}'] = int((a['best'].max(1) >= thr).sum())
        r[f'kgreedy{thr}'] = int((a['greedy'].max(1) >= thr).sum())
    if cfg['pa'] == 'IND':
        cyc_lock = 886.0 + a['ind_cyc']
        variants = {'main': cyc_lock}
    else:
        base = 886.0 + a['stage_cyc']
        ev = a['events']
        variants = {'main': base + ev * C_RS, 'copy0': base, 'copy4x_host5us': base + ev * 4 * C_RS}
    r['pop_cycles'] = {k: float(v.mean()) for k, v in variants.items()}
    k = r['k33337']; lo, hi = F.wilson(k, P)
    r['P_pop'] = k / P; r['P_pop_lo'] = lo; r['P_pop_hi'] = hi
    r['tts'] = {}
    for mx, tag in ((mux, ''), (1, '|R_engines')) if mux > 1 else ((mux, ''),):
        for vname, cyc in r['pop_cycles'].items():
            for fname, f in FREQS.items():
                host = (5e-6 * a['events']) if vname == 'copy4x_host5us' else 0.0
                t_pop = mx * cyc / f + host
                r['tts'][f'{vname}|{fname}{tag}'] = dict(
                    t_pop_ms=1e3 * t_pop, tts_ms=1e3 * F.tts_ms(t_pop, k / P) if k else math.inf,
                    tts_ms_lo=1e3 * F.tts_ms(t_pop, hi) if hi > 0 else math.inf,
                    tts_ms_hi=1e3 * F.tts_ms(t_pop, lo) if lo > 0 else math.inf, multiplex=mx)
    return r


def pa_score(r, thr=33337):
    lo, _ = F.wilson(r[f'k{thr}'], r['populations'])
    t_pop = r['pop_cycles']['main'] / FREQS['V80_250MHz']
    return 1e3 * F.tts_ms(t_pop, lo) if lo > 0 else math.inf


def pa_base():
    sel = load('a2_confirm')['selection']['onsager_overall']['cfg']
    return {k: sel[k] for k in ('family', 'T0', 'q', 'lam', 'ramp', 'T1')}


def stage_b():
    base = pa_base(); cs = pa_configs(base)
    log('stage b:', len(cs), 'PA configurations x 256 populations x 12', json.dumps(base), file_hashes())
    rows = []
    for i, cfg in enumerate(cs):
        t = time.time()
        a = run_pa(cfg, 256, 12, [BASE, 4, i], [BASE, 4, i, 1])
        r = pa_summary(cfg, a, 12); r['idx'] = i; r['sec'] = time.time() - t
        rows.append(r)
        log('B', i, key(cfg), f"Ppop337 {r['k33337']}/256 Ppop320 {r['k33320']} Ppop300 {r['k33300']} kg337 {r['kgreedy33337']} "
            f"single {r['single_replica_k33337']}/{r['single_replica_n']} copied/ev {r['mean_copied_per_event']:.2f} "
            f"t_pop {r['tts']['main|V80_250MHz']['t_pop_ms']:.3f} ms score {pa_score(r):.1f} ({r['sec']:.0f}s)")
        dump('b_pilot', dict(hashes=file_hashes(), base=base, rows=rows))
    pa_rows = [r for r in rows if r['cfg']['pa'] != 'IND']
    thr = 33337 if any(r['k33337'] > 0 for r in pa_rows) else 33320
    best = min(pa_rows, key=lambda r: (pa_score(r, thr), -r['mean_pop_max']))
    dump('b_pilot', dict(hashes=file_hashes(), base=base, rows=rows, selection=dict(cfg=best['cfg'], ranked_at=thr,
                                                                                   score_ms=pa_score(best, thr))))
    log('B SELECTION', json.dumps(best['cfg']), pa_score(best, thr))


def stage_hold_b(Rr=12, P_total=1024, chunk=256, stage=10, name='holdout_b'):
    sel = load('b_pilot')['selection']['cfg']
    ind = dict(sel, pa='IND', K=1, alpha=0.0)
    log(f'stage {name}: R={Rr}, {P_total} populations, PA', json.dumps(sel), file_hashes())
    out = {}
    for cfg in (sel, ind):
        t = time.time(); parts = []
        for ch in range(P_total // chunk):
            parts.append(run_pa(cfg, chunk, Rr, [BASE, stage, ch], [BASE, stage, ch, 1]))
        a = {k: (np.concatenate([p[k] for p in parts]) if isinstance(parts[0][k], np.ndarray) else parts[0][k])
             for k in parts[0]}
        r = pa_summary(cfg, a, Rr); r['sec'] = time.time() - t
        out[cfg['pa']] = r
        np.savez_compressed(RES / f'{name}_{cfg["pa"]}.npz', cut=a['cut'], best=a['best'], greedy=a['greedy'],
                            stage_cyc=a['stage_cyc'], ind_cyc=a['ind_cyc'])
        log(name, cfg['pa'], f"P_pop337 {r['k33337']}/{r['populations']} [{r['P_pop_lo']:.4f}, {r['P_pop_hi']:.4f}] "
            f"single {r['single_replica_k33337']}/{r['single_replica_n']} TTS V80 {r['tts']['main|V80_250MHz']['tts_ms']:.1f} ms "
            f"({r['sec']:.0f}s)")
        dump(name, dict(hashes=file_hashes(), selection=sel, results=out))


if __name__ == '__main__':
    st = sys.argv[1]
    {'a1': stage_a1, 'a2': stage_a2, 'c': stage_c, 'hold': stage_hold, 'b': stage_b,
     'hold_b': lambda: stage_hold_b(12, 1024, 256, 10, 'holdout_b'),
     'b48': lambda: stage_hold_b(48, 256, 64, 11, 'holdout_b48')}[st]()
