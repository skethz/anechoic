"""PROTOCOL_PUB.md runner: the K2000 algorithm comparison with every baseline at its published settings.
Usage: python3 run_pub.py run [workers] | summarize | selftest-small
Jobs (fresh seeds, never the frozen integer seeds):
  pilot  (SCA (STATICA), TEC on STATICA: 3 T_init per S; bSB, dSB: 5 Delta t per S): 256 runs,
         SeedSequence([20261008, 301, code, S index, k])
  final  every configuration of every baseline at every S: 256 runs, SeedSequence([20261008, 302, code, S index, k])
  ours   Onsager-kT and Onsager-online re-run at their frozen per-budget selections (results/S*.json of
         research/algorithm_compare_20261007), 256 runs, SeedSequence([20261008, 303, code, S index]) -- used only for the
         best-visited secondary estimator and as a reproducibility check; the primary values of ours are the frozen finals.
Sequential methods (SA (Neal), TEC (published)) are split into chunks of 32 runs; the chunks reproduce one 256-run call.
Raw results: results_pub/raw/*.json; summary: results_pub/summary.json and results_pub/S<S>.json (plot format)."""
import json
import math
import os
import sys
import time
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
import pub_methods as P  # noqa: E402
from common import M, ROOT  # noqa: E402

HERE = Path(__file__).resolve().parent
OUT = HERE / 'results_pub'
RAW = OUT / 'raw'
S_LIST = (250, 500, 1000, 2000, 4000)
B = 256
CHUNK = 32
SEED = 20261008
BASELINES = ('SA (Neal)', 'SCA (STATICA)', 'TEC (published)', 'TEC on STATICA', 'APC-SCA (published)', 'ReAIM ASA',
             'aSB', 'bSB', 'dSB')
PILOTED = ('SCA (STATICA)', 'TEC on STATICA', 'bSB', 'dSB')
OURS = ('TEC-T (ours)', 'Onsager SCA (ours)')
ORIG = ROOT / 'research/algorithm_compare_20261007/results'
if os.environ.get('PUB_TEST') == '1':     # smoke test only: other seed root, tiny budgets, separate output
    OUT = HERE / 'results_pub_TEST'; RAW = OUT / 'raw'; S_LIST = (16, 24); B = 64; SEED = 999


def jobs():
    out = []
    for si, S in enumerate(S_LIST):
        for name in BASELINES:
            code = P.CODE[name]
            for k in range(P.N_CFG[name]):
                cfg = P.cfg_of(name, S, k)
                kinds = ('pilot', 'final') if name in PILOTED else ('final',)
                for kind in kinds:
                    seed_words = [SEED, 301 if kind == 'pilot' else 302, code, si, k]
                    if cfg['family'] in P.SEQUENTIAL:
                        for c0 in range(0, B, CHUNK):
                            out.append(dict(name=name, S=S, si=si, k=k, kind=kind, cfg=cfg, seed=seed_words, chunk=[c0, c0 + CHUNK]))
                    else:
                        out.append(dict(name=name, S=S, si=si, k=k, kind=kind, cfg=cfg, seed=seed_words, chunk=None))
        if os.environ.get('PUB_TEST') != '1':
            for name in OURS:
                orig = json.loads((ORIG / f'S{S}.json').read_text())[name]['final']['cfg']
                out.append(dict(name=name, S=S, si=si, k=0, kind='ours', cfg=orig, seed=[SEED, 303, P.CODE[name], si], chunk=None))
    return out


def job_path(j):
    c = '' if j['chunk'] is None else f"_c{j['chunk'][0]:03d}"
    return RAW / f"{j['name'].replace(' ', '_').replace('(', '').replace(')', '')}_S{j['S']}_k{j['k']}_{j['kind']}{c}.json"


_J = None


def _init():
    global _J
    np.seterr(all='ignore')
    _J = M.m.load()


def _run(j):
    J, sumw = _J; t0 = time.time(); cfg = j['cfg']; S = j['S']; seed = np.random.SeedSequence(j['seed'])
    fam = cfg['family']
    if fam in P.SEQUENTIAL:
        s0, seeds = P.draw_sequential_inputs(seed, B, len(J))
        a, b = j['chunk']
        if fam == 'SA_neal':
            s, E = P.sa_neal_chunk(J, s0[a:b], seeds[a:b], S)
        else:
            s, E = P.tec_pub_chunk(J, s0[a:b], seeds[a:b], S)
        cuts = (sumw - E[-1]) / 2.0
        res = dict(E=np.rint(E).astype(np.int64).tolist(), cuts=cuts.tolist())   # H is integer-valued; exact
    elif j['kind'] == 'pilot':
        cuts, flips = P.run_pilot_sca(J, sumw, cfg, B, seed) if fam in ('plain', 'tec') else (None, None)
        if cuts is None:   # bSB / dSB pilots
            E, cuts, _ = P.run_vectorized(J, sumw, cfg, B, seed); flips = None
        res = dict(cuts=np.asarray(cuts).tolist(), mean_flips=None if flips is None else float(np.mean(flips)))
    else:
        E, cuts, _ = P.run_vectorized(J, sumw, cfg, B, seed)
        res = P.summarize_final(E, cuts, sumw, S)
    res.update(job={k: v for k, v in j.items()}, sec=time.time() - t0)
    path = job_path(j)
    path.write_text(json.dumps(res) + '\n')
    return str(path), res['sec']


def run(workers):
    RAW.mkdir(parents=True, exist_ok=True)
    todo = [j for j in jobs() if not job_path(j).exists()]
    todo.sort(key=lambda j: -j['S'] * (8 if j['cfg']['family'] in P.SEQUENTIAL else 1))   # long jobs first
    print(f'{len(todo)} jobs', flush=True)
    t0 = time.time()
    with ProcessPoolExecutor(max_workers=workers, initializer=_init) as ex:
        futs = [ex.submit(_run, j) for j in todo]
        for n, f in enumerate(as_completed(futs), 1):
            p, sec = f.result()
            print(f'[{n}/{len(todo)}] {time.time() - t0:7.0f}s {Path(p).name} ({sec:.0f}s)', flush=True)
    summarize()


# ------------------------------------------------------------------ selection rules (pilot data only) and summary
def statica_select(pilots, S):
    """STATICA's figure of merit: TTS = t ln(0.01)/ln(1-p) with STATICA's time model; all infinite -> highest mean cut."""
    sc = []
    for k, r in enumerate(pilots):
        cuts = np.array(r['cuts']); p = float((cuts >= 33000).mean())
        t = P.statica_t_ms(r['mean_flips'], S)
        sc.append((P.r99(p) * t if p > 0 else math.inf, -float(cuts.mean()), k))
    return min(sc)[2], [dict(k=k, tts_ms=s[0], mean_cut=-s[1]) for k, s in zip(range(len(sc)), sc)]


def sb_select(pil):
    """Goto et al. 2021: the best Delta t among five values for the problem, N_step optimised for the target:
    minimise over Delta t the minimum over S of S ln(0.01)/ln(1-p); ties -> highest pilot mean cut at that S."""
    rows = []
    for k in range(len(P.SB_DTS)):
        best = min(((P.mcs99(S, float((np.array(pil[(S, k)]['cuts']) >= 33000).mean())), S) for S in S_LIST))
        mc = float(np.mean(pil[(best[1], k)]['cuts']))
        rows.append(dict(k=k, dt=P.SB_DTS[k], min_mcs99=best[0], at_S=best[1], mean_cut_at_S=mc))
    sel = min(rows, key=lambda r: (r['min_mcs99'], -r['mean_cut_at_S'], r['k']))['k']
    return sel, rows


def load_raw():
    recs = {}
    for f in RAW.glob('*.json'):
        r = json.loads(f.read_text()); j = r['job']
        key = (j['name'], j['S'], j['k'], j['kind'])
        recs.setdefault(key, []).append(r)
    return recs


def merge_chunks(rs, S, sumw):
    rs = sorted(rs, key=lambda r: r['job']['chunk'][0])
    E = np.concatenate([np.array(r['E']) for r in rs], axis=1); cuts = np.concatenate([np.array(r['cuts']) for r in rs])
    assert E.shape[1] == B, (len(rs), E.shape)
    return P.summarize_final(E, cuts, sumw, S)


def summarize():
    J, sumw = M.m.load()
    recs = load_raw()
    finals, pilots = {}, {}
    for (name, S, k, kind), rs in recs.items():
        if kind == 'pilot':
            pilots[(name, S, k)] = rs[0]
        elif rs[0]['job']['chunk'] is not None:
            finals[(name, S, k, kind)] = merge_chunks(rs, S, sumw)
        else:
            finals[(name, S, k, kind)] = rs[0]
    sel = {}
    for S in S_LIST:
        for name in ('SCA (STATICA)', 'TEC on STATICA'):
            if all((name, S, k) in pilots for k in range(3)):
                sel[(name, S)] = statica_select([pilots[(name, S, k)] for k in range(3)], S)
    for name in ('bSB', 'dSB'):
        pil = {(S, k): pilots[(name, S, k)] for S in S_LIST for k in range(5) if (name, S, k) in pilots}
        if len(pil) == 5 * len(S_LIST):
            s_k, rows = sb_select(pil)
            for S in S_LIST:
                sel[(name, S)] = (s_k, rows)
    summary = dict(protocol='PROTOCOL_PUB.md', S_list=list(S_LIST), runs=B, methods={})
    plot = {S: {} for S in S_LIST}
    for name in BASELINES + OURS:
        summary['methods'][name] = {}
        for S in S_LIST:
            kind = 'ours' if name in OURS else 'final'
            ks = sorted(k for (n, s, k, kd) in finals if n == name and s == S and kd == kind)
            if not ks:
                continue
            main = sel[(name, S)][0] if (name, S) in sel else ks[0]
            entry = dict(selected_k=main, selection=None if (name, S) not in sel else sel[(name, S)][1],
                         configs={k: dict(cfg=finals[(name, S, k, kind)]['job']['cfg'] if 'job' in finals[(name, S, k, kind)] else
                                          P.cfg_of(name, S, k),
                                          **{x: finals[(name, S, k, kind)][x] for x in ('mean_cut', 'sd_cut', 'p33000',
                                             'p33000_wilson95', 'p33200', 'mcs99', 'best_visited', 'valid')}) for k in ks})
            summary['methods'][name][S] = entry
            if name not in OURS:
                f = dict(finals[(name, S, main, kind)]); f['cfg'] = entry['configs'][main]['cfg']
                plot[S][name] = dict(final=f, selected=main)
    if os.environ.get('PUB_TEST') != '1':
        for S in S_LIST:     # ours: the frozen finals (primary), unchanged
            orig = json.loads((ORIG / f'S{S}.json').read_text())
            for name in OURS:
                plot[S][name] = dict(final=orig[name]['final'], selected=orig[name]['selected'], source='frozen sweep (PROTOCOL.md)')
    OUT.mkdir(exist_ok=True)
    for S in S_LIST:
        (OUT / f'S{S}.json').write_text(json.dumps(plot[S]) + '\n')
    (OUT / 'summary.json').write_text(json.dumps(summary, indent=1, default=str) + '\n')
    table(summary, plot)


def table(summary, plot):
    lines = []
    for S in S_LIST:
        lines.append(f'=== S = {S}')
        for name, r in plot[S].items():
            f = r['final']
            bv = f.get('best_visited', {})
            lines.append(f"{name:22s} k={r.get('selected')} mean_cut={f['mean_cut']:9.1f} p33000={f['p33000']:.3f} "
                         f"mcs99={f['mcs99']:9.0f}  best-visited p={bv.get('p33000', float('nan')):.3f}  cfg={f.get('cfg')}")
    lines.append('=== steps to solution (min over S of S ln0.01/ln(1-p), final state)')
    for name in plot[S_LIST[0]]:
        v = [(plot[S][name]['final']['mcs99'], S) for S in S_LIST if name in plot[S]]
        lines.append(f'{name:22s} {min(v)[0]:9.0f} (at S = {min(v)[1]})')
    t = '\n'.join(lines)
    print(t)
    (OUT / 'table.txt').write_text(t + '\n')


def selftest_small():
    """Chunked sequential runs equal one unchunked call; job seeds are unique."""
    js = jobs(); seen = set()
    for j in js:
        key = (tuple(j['seed']), None if j['chunk'] is None else j['chunk'][0])
        assert key not in seen, key; seen.add(key)
    print(len(js), 'jobs, all seeds unique')


if __name__ == '__main__':
    cmd = sys.argv[1] if len(sys.argv) > 1 else 'selftest-small'
    if cmd == 'run':
        run(int(sys.argv[2]) if len(sys.argv) > 2 else 48)
    elif cmd == 'summarize':
        summarize()
    else:
        selftest_small()
