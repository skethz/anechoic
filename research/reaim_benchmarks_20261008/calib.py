"""EXPLORATORY pre-protocol calibration on synthetic instances only (never a G-set or TSPLIB benchmark instance).
Places the tuning grids of PROTOCOL.md. Stage 1: SA, plain SCA, ReAIM, aSB, bSB, dSB on broad log-spaced grids.
Stage 2: TEC, Onsager-kT, Onsager-online and APC-SCA around plain SCA's best (q, T0, T_fin) region.
Synthetic instances: GPP on random G(800, 19176) and planar-like 800-node graphs (research/gset_20261007/synth.py
generators, +1 weights); TSP on random Euclidean instances (n = 17, 29; integer distances, exact optimum by MILP).
Usage: python3 calib.py stage1|stage2|penalty [workers]. Results: calib/<stage>/<problem>/<method>.json"""
import itertools
import json
import math
import multiprocessing as mp
import os
import sys
import time
import zlib
from pathlib import Path

for _v in ('NUMBA_NUM_THREADS', 'OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS'):
    os.environ[_v] = '1'
sys.dont_write_bytecode = True
HERE = Path(__file__).resolve().parent
RES = HERE.parent
os.environ.setdefault('NUMBA_CACHE_DIR', str(HERE / '.numba_cache'))
for d in ('gset_20261007', 'ablation_20261005', 'theory_ideas_20261003', 'statica_reproduction_20261003'):
    sys.path.insert(0, str(RES / d))
sys.path.insert(0, str(HERE))
import numpy as np  # noqa: E402

import problems as PB  # noqa: E402
import solvers as SV  # noqa: E402

OUT = HERE / 'calib'
B_CAL = 16
S_CAL = {'gpp': 1024, 'tsp': 2048}
P_GPP = 4          # pre-declared penalty for the grids (stage 1/2); the 'penalty' stage checks P in {2, 4, 8}


# ------------------------------------------------------------------ synthetic instances
def syn_gpp(kind, seed, P=P_GPP):
    import synth
    g = synth.random_graph(800, 19176, False, seed) if kind == 'rand' else synth.planar_like(800, False, seed)
    return PB.gpp_from_edges(f'syn_gpp_{kind}{seed}_P{P}', g.N, g.ei, g.ej, g.w, P, ref=None, target=None)


def syn_tsp(n, seed):
    import tsplib
    rng = np.random.default_rng(seed)
    pts = rng.random((n, 2)) * 1000.0
    W = np.rint(np.hypot(*(pts[:, None, :] - pts[None, :, :]).transpose(2, 0, 1))).astype(np.int64)
    W[(W == 0) & ~np.eye(n, dtype=bool)] = 1
    L, _, _ = tsplib.exact_optimum(W)
    return PB.tsp_from_matrix(f'syn_tsp{n}_{seed}', W, float(W.max()), ref=L, target=int(math.floor(1.01 * L)))


INSTANCES = {'gpp': [('rand', 101), ('planar', 102)], 'tsp': [(17, 201), (29, 202)]}


def instance(prob, spec, P=P_GPP):
    return syn_gpp(*spec, P=P) if prob == 'gpp' else syn_tsp(*spec)


# ------------------------------------------------------------------ broad grids (field values in units of sigma_T)
def lg(a, b, k):
    return tuple(float(x) for x in np.geomspace(a, b, k))


def grid1(prob, method, P, S):
    a = P.sigma_T; st = P.stats(); N = P.N; out = []
    if method == 'SA':
        for t0, t1 in itertools.product(lg(0.125, 4, 6), lg(0.005, 0.32, 7)):
            if t1 < t0:
                out.append(dict(family='SA', T0=t0 * a, T1=t1 * a, S=S, rel=dict(T0=t0, T1=t1)))
    elif method == 'SCA':
        for q, t0, tf in itertools.product(lg(0.25, 16, 7), lg(0.25, 4, 5), (0.05, 0.15, 0.45)):
            if tf < t0:
                out.append(dict(family='plain', q=q * a, T0=t0 * a, tfin=tf * a, S=S, rel=dict(q=q, T0=t0, tfin=tf)))
    elif method == 'ReAIM ASA':
        ks_all = ((1, 2, 4, 8), (2, 4, 8, 16), (4, 8, 16, 32), (8, 16, 32, 64), (16, 32, 64, 128), (32, 64, 128, 256),
                  (64, 128, 256, 512))
        for ks, t1 in itertools.product(ks_all, (0.2, 0.1, 0.05, 0.025, 0.0125)):
            kk = tuple(int(min(N, k)) for k in ks)
            out.append(dict(family='ReAIM', kset=kk, T1=t1, S=S, rel=dict(kset=ks, T1=t1)))
    elif method in ('aSB', 'bSB', 'dSB'):
        dts = (0.125, 0.25, 0.5, 0.9, 1.25) if method == 'aSB' else (0.25, 0.5, 0.75, 1.0, 1.25)
        for d, x in itertools.product(dts, lg(0.25, 8, 6)):
            out.append(dict(family=method, dt=d, xi=x, S=S))
    return out


def grid1b(prob, method, P, S):
    """Extensions beyond stage-1 edges: SA T1 up; SB xi up; plain SCA q and T0 up."""
    a = P.sigma_T; N = P.N; out = []
    if method == 'SA':
        for t0, t1 in itertools.product(lg(0.5, 16, 6), (0.32, 0.64, 1.28)):
            if t1 < t0:
                out.append(dict(family='SA', T0=t0 * a, T1=t1 * a, S=S, rel=dict(T0=t0, T1=t1)))
    elif method == 'SCA':
        for q, t0, tf in itertools.product(lg(16, 128, 4), lg(0.5, 16, 6), (0.05, 0.15, 0.45, 1.35)):
            if tf < t0:
                out.append(dict(family='plain', q=q * a, T0=t0 * a, tfin=tf * a, S=S, rel=dict(q=q, T0=t0, tfin=tf)))
    elif method in ('aSB', 'bSB', 'dSB'):
        dts = (0.0625, 0.125, 0.25, 0.5) if method == 'aSB' else (0.125, 0.25, 0.5, 0.75, 1.0)
        for d, x in itertools.product(dts, lg(8, 64, 4)):
            out.append(dict(family=method, dt=d, xi=x, S=S))
    elif method == 'ReAIM ASA':
        for ks, t1 in itertools.product(((1, 1, 2, 2), (1, 2, 3, 4), (1, 2, 4, 8)), (0.4, 0.2, 0.1, 0.05, 0.025, 0.0125, 0.00625)):
            out.append(dict(family='ReAIM', kset=ks, T1=t1, S=S, rel=dict(kset=ks, T1=t1)))
    return out


def grid2b(prob, method, P, S):
    """Broad stage 2 for the corrected rules and APC-SCA: rule parameter x q x T0 (tfin per problem from stage 1)."""
    a = P.sigma_T; st = P.stats(); out = []
    tf = 0.15
    qs = lg(0.5, 64, 8); ts = lg(0.25, 8, 6)
    if method == 'TEC':
        for j, q, t in itertools.product((-4.0, -2.0, -1.0, -0.5, -0.25, 0.25), qs, ts):
            out.append(dict(family='tec', jv=j * a, q=q * a, T0=t * a, tfin=tf * a, S=S, rel=dict(jv=j, q=q, T0=t, tfin=tf)))
    elif method == 'Onsager-kT':
        for k, q, t in itertools.product(tuple(2.0 ** e for e in (-12, -10, -8, -6, -4, -2)), qs, ts):
            out.append(dict(family='tecT', kappa=k * st['kappa_factor'], q=q * a, T0=t * a, ramp=True, tfin=tf * a, S=S,
                            rel=dict(kappa=k, kappa_factor=st['kappa_factor'], q=q, T0=t, tfin=tf)))
    elif method == 'Onsager-online':
        for lam, q, t in itertools.product(tuple(2.0 ** e for e in (-12, -10, -8, -6, -4, -2)), qs, ts):
            out.append(dict(family='onsager', lam=lam * st['lam_factor'], q=q * a, T0=t * a, ramp=True, tfin=tf * a, S=S,
                            rel=dict(lam=lam, lam_factor=st['lam_factor'], q=q, T0=t, tfin=tf)))
    elif method == 'APC-SCA':
        for x, r, y, t in itertools.product((2.0, 8.0, 32.0, 128.0), (0.9, 0.97), (0.25, 1.0, 4.0, 16.0), ts):
            if y < x:
                out.append(dict(family='apc', q_reset=x * a, r_q=r, q_lim=y * a, T0=t * a, tfin=tf * a, S=S,
                                rel=dict(q_reset=x, r_q=r, q_lim=y, T0=t, tfin=tf)))
    return out


def grid2c(prob, method, P, S):
    """Extension beyond stage-2b edges (larger q, T0, more negative J_v, larger q_reset / r_q) for all SCA-family rules,
    and stage-1 edges (SA T1, SB dt) not yet covered."""
    a = P.sigma_T; st = P.stats(); out = []
    qs = lg(8, 128, 5); ts = lg(2, 32, 5); tfs = (0.15, 0.6)
    if method == 'SCA':
        for q, t, tf in itertools.product(qs, ts, tfs):
            out.append(dict(family='plain', q=q * a, T0=t * a, tfin=tf * a, S=S, rel=dict(q=q, T0=t, tfin=tf)))
    elif method == 'TEC':
        for j, q, t, tf in itertools.product((-16.0, -8.0, -4.0), qs, ts, tfs):
            out.append(dict(family='tec', jv=j * a, q=q * a, T0=t * a, tfin=tf * a, S=S, rel=dict(jv=j, q=q, T0=t, tfin=tf)))
    elif method == 'Onsager-kT':
        for k, q, t, tf in itertools.product(tuple(2.0 ** e for e in (-8, -6, -4, -2, 0)), qs, ts, tfs):
            out.append(dict(family='tecT', kappa=k * st['kappa_factor'], q=q * a, T0=t * a, ramp=True, tfin=tf * a, S=S,
                            rel=dict(kappa=k, kappa_factor=st['kappa_factor'], q=q, T0=t, tfin=tf)))
    elif method == 'Onsager-online':
        for lam, q, t, tf in itertools.product(tuple(2.0 ** e for e in (-8, -6, -4, -2, 0)), qs, ts, tfs):
            out.append(dict(family='onsager', lam=lam * st['lam_factor'], q=q * a, T0=t * a, ramp=True, tfin=tf * a, S=S,
                            rel=dict(lam=lam, lam_factor=st['lam_factor'], q=q, T0=t, tfin=tf)))
    elif method == 'APC-SCA':
        for x, r, y, t in itertools.product((64.0, 128.0, 256.0, 512.0), (0.97, 0.99), (0.125, 0.25, 1.0, 4.0, 8.0), lg(0.5, 8, 5)):
            out.append(dict(family='apc', q_reset=x * a, r_q=r, q_lim=y * a, T0=t * a, tfin=0.15 * a, S=S,
                            rel=dict(q_reset=x, r_q=r, q_lim=y, T0=t, tfin=0.15)))
    elif method == 'SA':
        for t0, t1 in itertools.product(lg(1, 8, 4), (0.64, 1.28, 2.56)):
            if t1 < t0:
                out.append(dict(family='SA', T0=t0 * a, T1=t1 * a, S=S, rel=dict(T0=t0, T1=t1)))
    elif method in ('bSB', 'dSB', 'aSB'):
        for d, x in itertools.product((0.03125, 0.0625, 0.125), lg(2, 32, 5)):
            out.append(dict(family=method, dt=d, xi=x, S=S))
    return out


def grid2(prob, method, P, S, base):
    """Stage 2 around plain SCA's best (q, T0, tfin) (relative units): q x {1/2, 1, 2}, T0 x {1/2, 1, 2}, tfin fixed."""
    a = P.sigma_T; st = P.stats(); out = []
    qb, tb, tf = base['q'], base['T0'], base['tfin']
    qs = (qb / 2, qb, qb * 2); ts = tuple(t for t in (tb / 2, tb, tb * 2) if t > tf)
    if method == 'TEC':
        for j, q, t in itertools.product((-1.0, -0.5, -0.25, -0.12, -0.06, 0.06, 0.12, 0.25), qs, ts):
            out.append(dict(family='tec', jv=j * a, q=q * a, T0=t * a, tfin=tf * a, S=S, rel=dict(jv=j, q=q, T0=t, tfin=tf)))
    elif method == 'Onsager-kT':
        for k, q, t in itertools.product(lg(2 ** -12, 2 ** 1, 14)[::2], qs, ts):
            out.append(dict(family='tecT', kappa=k * st['kappa_factor'], q=q * a, T0=t * a, ramp=True, tfin=tf * a, S=S,
                            rel=dict(kappa=k, kappa_factor=st['kappa_factor'], q=q, T0=t, tfin=tf)))
    elif method == 'Onsager-online':
        for lam, q, t in itertools.product(lg(2 ** -12, 2 ** 1, 14)[::2], qs, ts):
            out.append(dict(family='onsager', lam=lam * st['lam_factor'], q=q * a, T0=t * a, ramp=True, tfin=tf * a, S=S,
                            rel=dict(lam=lam, lam_factor=st['lam_factor'], q=q, T0=t, tfin=tf)))
    elif method == 'APC-SCA':
        for x, r, y, t in itertools.product((qb * 2, qb * 4, qb * 8), (0.9, 0.97), (qb / 4, qb / 2, qb), ts):
            out.append(dict(family='apc', q_reset=x * a, r_q=r, q_lim=y * a, T0=t * a, tfin=tf * a, S=S,
                            rel=dict(q_reset=x, r_q=r, q_lim=y, T0=t, tfin=tf)))
    return out


# ------------------------------------------------------------------ jobs
_CACHE = {}


def get_problem(prob, spec, P=P_GPP):
    key = (prob, spec, P)
    if key not in _CACHE:
        _CACHE[key] = instance(prob, spec, P)
    return _CACHE[key]


def run_point(P, cfg, seed):
    s, flips, finite = SV.run_method(P, cfg, B_CAL, seed)
    if P.kind == 'gpp':
        si = s[:, P.meta['ei']]; sj = s[:, P.meta['ej']]
        cut = ((si != sj).astype(np.int64) * P.meta['w'][None, :]).sum(1)
        M = np.rint(s.astype(np.float64).sum(1)).astype(np.int64)
        ok = np.ones(len(s), bool) if finite is None else np.asarray(finite, bool)
        return dict(cut=cut.tolist(), M=M.tolist(), ok=ok.tolist(), flips=None if flips is None else float(np.mean(flips)))
    ev = P.evaluate(s)
    ok = np.ones(len(s), bool) if finite is None else np.asarray(finite, bool)
    return dict(L=ev['value'].tolist(), feasible=(ev['feasible'] & ok).tolist(), rows_ok=ev['rows_ok'].tolist(),
                cols_ok=ev['cols_ok'].tolist(), flips=None if flips is None else float(np.mean(flips)))


def job(args):
    stage, prob, spec, method, cfgs, P_pen = args
    path = OUT / stage / prob / f"{prob}_{'_'.join(map(str, spec))}_P{P_pen}" / f"{method.replace(' ', '_')}.json"
    if path.exists():
        return None
    t0 = time.time()
    P = get_problem(prob, spec, P_pen)
    recs = []
    for i, cfg in enumerate(cfgs):
        r = run_point(P, cfg, np.random.SeedSequence([777, zlib.crc32(f'{stage}|{prob}|{spec}|{method}'.encode()), i]))
        recs.append(dict(cfg=cfg, **r))
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(dict(stage=stage, problem=prob, instance=P.name, spec=spec, method=method, ref=P.meta.get('ref'),
                                    stats=P.stats(), points=recs, sec=time.time() - t0)) + '\n')
    return f'{stage} {P.name} {method}: {len(cfgs)} points {time.time() - t0:.0f}s'


def chunks(lst, k):
    return [lst[i::k] for i in range(k)]


def main():
    stage = sys.argv[1]; W = int(sys.argv[2]) if len(sys.argv) > 2 else 32
    jobs = []
    if stage == 'stage1':
        for prob, specs in INSTANCES.items():
            for spec in specs:
                P = get_problem(prob, spec)
                for method in ('SA', 'SCA', 'ReAIM ASA', 'aSB', 'bSB', 'dSB'):
                    cf = grid1(prob, method, P, S_CAL[prob])
                    nchunk = max(1, len(cf) // 12)
                    for c, part in enumerate(chunks(cf, nchunk)):
                        jobs.append(('stage1', prob, spec, f'{method}#{c}', part, P_GPP))
    elif stage == 'stage1b':
        for prob, specs in INSTANCES.items():
            for spec in specs:
                P = get_problem(prob, spec)
                for method in ('SA', 'SCA', 'ReAIM ASA', 'aSB', 'bSB', 'dSB'):
                    cf = grid1b(prob, method, P, S_CAL[prob])
                    nchunk = max(1, len(cf) // 12)
                    for c, part in enumerate(chunks(cf, nchunk)):
                        jobs.append(('stage1b', prob, spec, f'{method}#{c}', part, P_GPP))
    elif stage == 'stage2b':
        for prob, specs in INSTANCES.items():
            for spec in specs:
                P = get_problem(prob, spec)
                for method in ('TEC', 'Onsager-kT', 'Onsager-online', 'APC-SCA'):
                    cf = grid2b(prob, method, P, S_CAL[prob])
                    nchunk = max(1, len(cf) // 12)
                    for c, part in enumerate(chunks(cf, nchunk)):
                        jobs.append(('stage2b', prob, spec, f'{method}#{c}', part, P_GPP))
    elif stage == 'stage2c':
        for prob, specs in INSTANCES.items():
            for spec in specs:
                P = get_problem(prob, spec)
                for method in ('SCA', 'TEC', 'Onsager-kT', 'Onsager-online', 'APC-SCA', 'SA', 'bSB', 'dSB', 'aSB'):
                    cf = grid2c(prob, method, P, S_CAL[prob])
                    nchunk = max(1, len(cf) // 12)
                    for c, part in enumerate(chunks(cf, nchunk)):
                        jobs.append(('stage2c', prob, spec, f'{method}#{c}', part, P_GPP))
    elif stage == 'stage2':
        base = json.loads((OUT / 'stage2_base.json').read_text())
        for prob, specs in INSTANCES.items():
            for spec in specs:
                P = get_problem(prob, spec)
                for method in ('TEC', 'Onsager-kT', 'Onsager-online', 'APC-SCA'):
                    cf = grid2(prob, method, P, S_CAL[prob], base[prob])
                    nchunk = max(1, len(cf) // 12)
                    for c, part in enumerate(chunks(cf, nchunk)):
                        jobs.append(('stage2', prob, spec, f'{method}#{c}', part, P_GPP))
    elif stage == 'penalty':
        base = json.loads((OUT / 'stage2_base.json').read_text())
        for spec in INSTANCES['gpp']:
            for Ppen in (2, 4, 8):
                P = get_problem('gpp', spec, Ppen)
                bs = base['gpp_best_cfgs']
                for method, rel in bs.items():
                    cfg = dict(rel['cfg'])
                    jobs.append(('penalty', 'gpp', spec, method, [cfg], Ppen))
    print(f'{stage}: {len(jobs)} jobs, {W} workers', flush=True)
    t0 = time.time()
    with mp.Pool(W) as pool:
        for k, msg in enumerate(pool.imap_unordered(job, jobs, chunksize=1)):
            if msg:
                print(f'[{k + 1}/{len(jobs)} {time.time() - t0:.0f}s] {msg}', flush=True)
    print(f'{stage} done {time.time() - t0:.0f}s', flush=True)


if __name__ == '__main__':
    main()
