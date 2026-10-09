"""PROTOCOL.md runner. Usage: python3 run_gset.py A | C | B
  A: per-budget pilot (64/grid point) -> selection by pilot mean cut -> 256-run final, four engine rules, all 51 instances
  C: TTS confirmations (1024 paired runs at S*_1 and S*_12 per engine rule and instance), after A
  B: per-budget pilot/final for SA, APC-SCA, ReAIM ASA, aSB, bSB, dSB on the 12-instance subset
Results: results/budget/G<g>/<method>_S<S>.json, results/conf/G<g>/<method>_<E1|E12>.json (resumable; existing files skipped).
6 single-threaded worker processes."""
import os

for _v in ('NUMBA_NUM_THREADS', 'OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'VECLIB_MAXIMUM_THREADS', 'MKL_NUM_THREADS'):
    os.environ[_v] = '1'

import itertools  # noqa: E402
import json  # noqa: E402
import math  # noqa: E402
import multiprocessing as mp  # noqa: E402
import sys  # noqa: E402
import time  # noqa: E402

import numpy as np  # noqa: E402

import engine as E  # noqa: E402
import others as O  # noqa: E402

ROOT = E.ROOT
OUT = ROOT / 'results'
S_LIST = (250, 500, 1000, 2000, 4000)
B_PILOT, B_FINAL, B_CONF = 64, 256, 1024
SEED = 20261007
WORKERS = 6
SUBSET = (1, 6, 11, 14, 18, 22, 27, 32, 35, 39, 43, 51)
ALL = tuple(range(1, 48)) + (51, 52, 53, 54)
METHODS = ('SA', 'SCA', 'TEC', 'APC-SCA', 'ReAIM ASA', 'aSB', 'bSB', 'dSB', 'Onsager-kT', 'Onsager-online')
ENGINE_RULES = ('SCA', 'TEC', 'Onsager-kT', 'Onsager-online')
OTHERS = ('SA', 'APC-SCA', 'ReAIM ASA', 'aSB', 'bSB', 'dSB')
OBJECTIVES = ('E1', 'E12')
ENGINE_FAMILIES = ('plain', 'tec', 'tecT', 'onsager', 'apc')

P_PLUS = set(range(1, 6)) | set(range(22, 27)) | set(range(43, 48))
Q_PLUS = set(range(14, 18)) | set(range(35, 39)) | set(range(51, 55))
SHARED = {
    'P+': dict(plain_q=(1.08, 1.44, 2.16, 2.88), plain_T0=(0.9, 1.35, 1.8, 2.7), q=(1.44, 2.16), T0=(1.35, 1.8), tfin=0.672,
               jv=(-0.36, -0.18, -0.09, 0.09), q_reset=(2.88, 4.32), q_lim=(1.44, 2.16)),
    'Q+': dict(plain_q=(0.72, 1.08, 1.44, 2.16), plain_T0=(0.45, 0.9, 1.35, 1.8), q=(1.08, 1.44), T0=(0.9, 1.35), tfin=0.672,
               jv=(-0.36, -0.18, -0.09, 0.09), q_reset=(2.16, 2.88), q_lim=(1.08, 1.44)),
    'M+-': dict(plain_q=(0.045, 0.09, 0.18, 0.36), plain_T0=(0.45, 0.67, 0.9, 1.35), q=(0.09, 0.18), T0=(0.45, 0.9), tfin=0.112,
                jv=(-0.09, -0.045, 0.045, 0.09), q_reset=(0.18, 0.36), q_lim=(0.045, 0.09)),
}
KAPPA = (0.125, 0.25, 0.5, 1.0)
LAMBDA = (0.125, 0.25, 0.5, 1.0)
SA_T0 = (0.447, 0.671, 1.118, 1.789)
SA_T1 = (0.0224, 0.0447, 0.0671, 0.1118)
REAIM_K = ((32, 64, 128, 256), (64, 128, 256, 512), (128, 256, 512, 1024), (256, 512, 1024, 2000))
REAIM_T1 = (0.2, 0.1, 0.05, 0.025)

BKV = json.loads((ROOT / 'bkv.json').read_text())['instances']


def wclass(g):
    return 'P+' if g in P_PLUS else 'Q+' if g in Q_PLUS else 'M+-'


def grid(method, G, g, S):
    """The 16-point grid of PROTOCOL.md for one method on instance g (Graph G), absolute values, with 'rel' recorded."""
    c = SHARED[wclass(g)]; a = G.sigma; tf = c['tfin'] * a; P = itertools.product
    out = []
    if method == 'SCA':
        for q, t in P(c['plain_q'], c['plain_T0']):
            out.append(dict(family='plain', q=q * a, T0=t * a, tfin=tf, S=S, rel=dict(q=q, T0=t, tfin=c['tfin'])))
    elif method == 'TEC':
        for j, q, t in P(c['jv'], c['q'], c['T0']):
            out.append(dict(family='tec', jv=j * a, q=q * a, T0=t * a, tfin=tf, S=S, rel=dict(jv=j, q=q, T0=t, tfin=c['tfin'])))
    elif method == 'Onsager-kT':
        for k, q, t in P(KAPPA, c['q'], c['T0']):
            out.append(dict(family='tecT', kappa=k, q=q * a, T0=t * a, ramp=True, tfin=tf, S=S,
                            rel=dict(kappa=k, q=q, T0=t, tfin=c['tfin'])))
    elif method == 'Onsager-online':
        for lam, q, t in P(LAMBDA, c['q'], c['T0']):
            out.append(dict(family='onsager', lam=lam * G.lam_factor, q=q * a, T0=t * a, ramp=True, tfin=tf, S=S,
                            rel=dict(lam_k2000=lam, lam_factor=G.lam_factor, q=q, T0=t, tfin=c['tfin'])))
    elif method == 'APC-SCA':
        for x, r, y, t in P(c['q_reset'], (0.9, 0.97), c['q_lim'], c['T0']):
            out.append(dict(family='apc', q_reset=x * a, r_q=r, q_lim=y * a, T0=t * a, tfin=tf, S=S,
                            rel=dict(q_reset=x, r_q=r, q_lim=y, T0=t, tfin=c['tfin'])))
    elif method == 'SA':
        for t0, t1 in P(SA_T0, SA_T1):
            out.append(dict(family='SA', T0=t0 * a, T1=t1 * a, S=S, rel=dict(T0=t0, T1=t1)))
    elif method == 'ReAIM ASA':
        for ks, t1 in P(REAIM_K, REAIM_T1):
            kk = tuple(int(min(G.N, max(1, round(k * G.N / 2000)))) for k in ks)
            out.append(dict(family='ReAIM', kset=kk, T1=t1, S=S, rel=dict(kset_n2000=ks, T1=t1)))
    elif method == 'aSB':
        for d, x in P((0.25, 0.5, 0.9, 1.25), (0.5, 0.75, 1.0, 1.5)):
            out.append(dict(family='aSB', dt=d, xi=x, S=S))
    elif method in ('bSB', 'dSB'):
        for d, x in P((0.5, 0.75, 1.0, 1.25), (0.5, 0.75, 1.0, 2.0)):
            out.append(dict(family=method, dt=d, xi=x, S=S))
    assert len(out) == 16, (method, len(out))
    return out


_GRAPHS = {}


def graph(g):
    if g not in _GRAPHS:
        _GRAPHS[g] = E.Graph(f'G{g}')
    return _GRAPHS[g]


def run_method(G, cfg, B, seed):
    """Returns exact integer cuts (B,), flips per run (engine families) or None, finite flags or None."""
    rng = np.random.default_rng(seed)
    fam = cfg['family']; flips = None; finite = None
    if fam in ENGINE_FAMILIES:
        s0 = rng.choice(np.array([-1.0, 1.0], dtype=np.float32), size=(B, G.N))
        s, flips = E.run(G, s0, cfg, rng, cfg['tfin'])
    elif fam == 'SA':
        s = O.sa(G, B, cfg, rng)
    elif fam in ('bSB', 'dSB'):
        s = O.sb(G, B, cfg, rng)
    elif fam == 'aSB':
        s, finite = O.asb(G, B, cfg, rng)
    elif fam == 'ReAIM':
        s = O.reaim(G, B, cfg, rng)
    else:
        raise ValueError(fam)
    return G.cut(s), flips, finite


def wilson(k, n, z=1.959963984540054):
    if n == 0:
        return (0.0, 1.0)
    p = k / n; d = 1 + z * z / n
    c = (p + z * z / (2 * n)) / d; h = z * math.sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / d
    return (max(0.0, c - h), min(1.0, c + h))


def mcs99(S, p):
    return math.inf if p <= 0 else (float(S) if p >= 1 else S * math.log(0.01) / math.log1p(-p))


def stats(cuts, flips, finite, g, S):
    b = BKV[f'G{g}']; tgt = b['target']
    ok = np.ones(len(cuts), bool) if finite is None else np.asarray(finite, bool)
    succ = (cuts >= tgt) & ok
    k = int(succ.sum()); n = len(cuts); p = k / n
    cv = cuts[ok]
    out = dict(runs=n, n_invalid=int((~ok).sum()), target=tgt, BKV=b['BKV'], k_target=k, p_target=p,
               p_target_wilson95=wilson(k, n), k_bkv=int(((cuts >= b['BKV']) & ok).sum()),
               mean_cut=float(cv.mean()) if len(cv) else None, sd_cut=float(cv.std(ddof=1)) if len(cv) > 1 else None,
               max_cut=int(cv.max()) if len(cv) else None, mcs99=mcs99(S, p))
    if g == 23:
        out['k_target_13221'] = int(((cuts >= 13221) & ok).sum())
    if flips is not None:
        out['mean_flips'] = float(np.mean(flips))
        tt = E.t_trial_ms(np.asarray(flips), S)
        out['t_trial_ms'] = float(E.t_trial_ms(float(np.mean(flips)), S))
        out['tts1_ms'] = out['t_trial_ms'] * E.r99(p)
        tr = E.expected_max(tt, E.E_ENGINES); P = 1 - (1 - p) ** E.E_ENGINES
        out['t_round12_ms'] = tr; out['P_round12'] = P
        out['tts12_ms'] = tr if P >= 0.995 else tr * E.r99(P)
    return out


def budget_path(g, method, S):
    return OUT / 'budget' / f'G{g}' / f"{method.replace(' ', '_')}_S{S}.json"


def job_budget(args):
    g, method, S = args
    path = budget_path(g, method, S)
    if path.exists():
        return None
    t0 = time.time()
    G = graph(g); mi = METHODS.index(method); si = S_LIST.index(S)
    gr = grid(method, G, g, S)
    tgt = BKV[f'G{g}']['target']
    pilot = []
    for gi, cfg in enumerate(gr):
        cuts, flips, finite = run_method(G, cfg, B_PILOT, np.random.SeedSequence([SEED, 1, g, mi, si, gi]))
        valid = finite is None or bool(np.all(finite))
        pilot.append(dict(cfg=cfg, valid=valid, mean_cut=float(cuts.mean()) if valid else -math.inf,
                          p_target=float((cuts >= tgt).mean()) if valid else 0.0,
                          mean_flips=float(np.mean(flips)) if flips is not None else None))
    ok = [i for i, r in enumerate(pilot) if r['valid']]
    rec = dict(instance=f'G{g}', method=method, method_index=mi, S=S, wclass=wclass(g), sigma=G.sigma, N=G.N, m=G.m,
               lam_factor=G.lam_factor, pilot=pilot)
    if ok:
        sel = max(ok, key=lambda i: (pilot[i]['mean_cut'], -i))
        cuts, flips, finite = run_method(G, gr[sel], B_FINAL, np.random.SeedSequence([SEED, 2, g, mi, si]))
        rec.update(selected=sel, selected_cfg=gr[sel], final=stats(cuts, flips, finite, g, S), cuts=cuts.tolist(),
                   flips=None if flips is None else np.asarray(flips).tolist(),
                   finite=None if finite is None else np.asarray(finite).tolist())
    else:
        rec.update(selected=None, selected_cfg=None, final=None)
    rec['sec'] = time.time() - t0
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_suffix('.tmp'); tmp.write_text(json.dumps(rec) + '\n'); tmp.rename(path)
    f = rec['final']
    return (f"G{g:<3d} {method:15s} S={S:<5d} sel={rec['selected']} pilot={pilot[rec['selected']]['mean_cut']:.1f} | "
            f"final mean={f['mean_cut']:.1f} max={f['max_cut']} p={f['p_target']:.3f} bkv={f['k_bkv']} mcs99={f['mcs99']:.0f}"
            + (f" tts1={f['tts1_ms']:.4g}ms tts12={f['tts12_ms']:.4g}ms" if 'tts1_ms' in f else '') + f" ({rec['sec']:.0f}s)"
            if f else f"G{g} {method} S={S}: no valid configuration ({rec['sec']:.0f}s)")


def conf_path(g, method, obj):
    return OUT / 'conf' / f'G{g}' / f"{method.replace(' ', '_')}_{obj}.json"


def choose_S(g, method, obj):
    """argmin over S of the objective's TTS on the 256-run finals (ties -> smaller S). None if infinite everywhere."""
    best = None
    for S in S_LIST:
        d = json.loads(budget_path(g, method, S).read_text())
        if d['final'] is None:
            continue
        v = d['final']['tts1_ms' if obj == 'E1' else 'tts12_ms']
        if math.isfinite(v) and (best is None or v < best[0]):
            best = (v, S, d['selected_cfg'])
    return best


def job_conf(args):
    g, method, obj = args
    path = conf_path(g, method, obj)
    if path.exists():
        return None
    t0 = time.time()
    ch = choose_S(g, method, obj)
    rec = dict(instance=f'G{g}', method=method, objective=obj)
    if ch is None:
        rec.update(S=None, cfg=None, final=None, note='TTS infinite at every budget in the 256-run finals: not confirmed')
    else:
        _, S, cfg = ch
        G = graph(g)
        cuts, flips, finite = run_method(G, cfg, B_CONF, np.random.SeedSequence([SEED, 3, g, OBJECTIVES.index(obj)]))
        rec.update(S=S, cfg=cfg, tts_selection_256=ch[0], final=stats(cuts, flips, finite, g, S), cuts=cuts.tolist(),
                   flips=np.asarray(flips).tolist())
    rec['sec'] = time.time() - t0
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_suffix('.tmp'); tmp.write_text(json.dumps(rec) + '\n'); tmp.rename(path)
    f = rec['final']
    return (f"CONF G{g:<3d} {method:15s} {obj:3s} S={rec['S']} p={f['p_target']:.3f} bkv={f['k_bkv']} "
            f"tts1={f['tts1_ms']:.4g}ms tts12={f['tts12_ms']:.4g}ms ({rec['sec']:.0f}s)" if f else
            f"CONF G{g} {method} {obj}: not confirmed (infinite)")


def order(jobs, cost):
    return sorted(jobs, key=lambda j: -cost(j))


def main():
    part = sys.argv[1]
    if part == 'A':
        sub = order([(g, m, S) for g in SUBSET for m in ENGINE_RULES for S in S_LIST], lambda j: graph_n(j[0]) * j[2])
        rest = order([(g, m, S) for g in ALL if g not in SUBSET for m in ENGINE_RULES for S in S_LIST],
                     lambda j: graph_n(j[0]) * j[2])
        jobs, fn = sub + rest, job_budget
    elif part == 'C':
        jobs = order([(g, m, o) for g in SUBSET + tuple(x for x in ALL if x not in SUBSET) for m in ENGINE_RULES
                      for o in OBJECTIVES], lambda j: 0)
        jobs, fn = jobs, job_conf
    elif part == 'B':
        jobs = order([(g, m, S) for g in SUBSET for m in OTHERS for S in S_LIST], lambda j: graph_n(j[0]) * j[2])
        fn = job_budget
    else:
        raise SystemExit(__doc__)
    print(f'part {part}: {len(jobs)} jobs, {WORKERS} workers, start {time.strftime("%Y-%m-%d %H:%M:%S")}', flush=True)
    t0 = time.time(); done = 0
    with mp.Pool(WORKERS) as pool:
        for msg in pool.imap_unordered(fn, jobs, chunksize=1):
            done += 1
            if msg:
                print(f'[{done}/{len(jobs)} {time.time() - t0:.0f}s] {msg}', flush=True)
    print(f'part {part} complete {time.strftime("%Y-%m-%d %H:%M:%S")} ({time.time() - t0:.0f}s)', flush=True)


def graph_n(g):
    return 800 if g <= 21 else 2000 if g <= 42 else 1000


if __name__ == '__main__':
    main()
