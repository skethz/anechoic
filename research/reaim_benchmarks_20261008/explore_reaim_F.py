"""EXPLORATORY, POST HOC (not part of the frozen protocol): ReAIM ASA with F = min (greedy-leaning q from the |N| FIFO).
ReAIM's Step 4 (Sec. IV-B) selects F per problem between max (SA-like) and min (greedy-like) and reports that TSP leans
towards min, while the implementation inherited from the K2000 / G-set studies (others.py, solvers.py) fixes F = max.
Same 32-point grid, pilot 64 runs per point, selection by pilot mean quality, 256-run final, at ReAIM's budgets
(GPP S = 4096, TSP S = 8192), fresh seeds SeedSequence([20261008, 9, p, i, grid index]) / [.., 10, p, i].
Also reports "Step-4 selection": per instance, the F whose pilot-selected configuration has the higher pilot mean
quality (protocol F = max pilot vs this F = min pilot). Results: explore/reaim_F/<problem>/<instance>.json"""
import json
import math
import multiprocessing as mp
import os
import sys
import time
from pathlib import Path

for _v in ('NUMBA_NUM_THREADS', 'OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS'):
    os.environ[_v] = '1'
sys.dont_write_bytecode = True
HERE = Path(__file__).resolve().parent
os.environ.setdefault('NUMBA_CACHE_DIR', str(HERE / '.numba_cache'))
sys.path.insert(0, str(HERE))
import numpy as np  # noqa: E402

import grids  # noqa: E402
import run_bench as RB  # noqa: E402
import solvers as SV  # noqa: E402

OUT = HERE / 'explore' / 'reaim_F'
S_OF = {'gpp': 4096, 'tsp': 8192}


def reaim_F(P, B, cfg, rng, F):
    """solvers.reaim with the |N|-FIFO reduction F (np.max in the protocol)."""
    S = cfg['S']; kset = cfg['kset']; t0, t1 = 1.0, cfg['T1']; it_trial, it_run = 32, 96
    N = P.N; R = len(kset)
    alpha = (t1 / t0) ** (1.0 / (S - 1))
    x = rng.choice(np.array([-1.0, 1.0], dtype=np.float32), size=(B, N))
    cur = SV.Batch(P, x, np.full((B, 20), N, dtype=np.int64), F)
    t, T = 0, t0
    while t < S:
        n = min(it_trial, S - t)
        rep = SV.Batch(P, np.repeat(cur.x, R, 0), np.repeat(cur.fifo, R, 0), F, h=np.repeat(cur.h, R, 0),
                       M=np.repeat(cur.M, R, 0))
        kvec = np.tile(np.array(kset), B)
        for u in range(n):
            rep.step(kvec, T, rng); T *= alpha
        pick = rep.energy().reshape(B, R).argmin(1)
        t += n
        sel = np.arange(B) * R + pick
        cur = SV.Batch(P, rep.x[sel].copy(), rep.fifo[sel].copy(), F, h=rep.h[sel].copy(), M=rep.M[sel].copy())
        kbest = np.array(kset)[pick]
        n = min(it_run, S - t)
        for u in range(n):
            cur.step(kbest, T, rng); T *= alpha
        t += n
    return cur.x


def job(args):
    prob, inst = args
    tag = f'G{inst}' if prob == 'gpp' else inst
    path = OUT / prob / f'{tag}.json'
    if path.exists():
        return None
    t0 = time.time()
    P = RB.problem(prob, inst); S = S_OF[prob]; pc = RB.PCODE[prob]; ic = RB.icode(prob, inst)
    gr = grids.grid(prob, 'ReAIM ASA', P, S)
    pilot = []
    for gi, cfg in enumerate(gr):
        rng = np.random.default_rng(np.random.SeedSequence([20261008, 9, pc, ic, gi]))
        s = reaim_F(P, 64, cfg, rng, np.min)
        ev = P.evaluate(s)
        pilot.append(dict(mean_quality=float(ev['quality'].mean()), p_feasible=float(ev['feasible'].mean())))
    sel = max(range(len(gr)), key=lambda i: (pilot[i]['mean_quality'], -i))
    rng = np.random.default_rng(np.random.SeedSequence([20261008, 10, pc, ic]))
    s = reaim_F(P, 256, gr[sel], rng, np.min)
    ev = P.evaluate(s)
    prot = json.loads((HERE / 'results' / prob / tag / f'ReAIM_ASA_S{S}.json').read_text())
    prot_pilot = prot['pilot'][prot['selected']]['mean_quality']
    rec = dict(problem=prob, instance=P.name, S=S, F='min', selected=sel, selected_cfg=gr[sel], pilot=pilot,
               mean_quality=float(ev['quality'].mean()), p_feasible=float(ev['feasible'].mean()),
               p_target=float(ev['success'].mean()), k_opt=int(ev['opt'].sum()),
               best=int(ev['value'][ev['feasible']].min()) if ev['feasible'].any() else None,
               protocol_F_max=dict(pilot_best=prot_pilot, final_mean_quality=prot['final']['mean_quality'],
                                   final_p_feasible=prot['final']['p_feasible'], final_p_target=prot['final']['p_target']),
               step4_choice='min' if pilot[sel]['mean_quality'] > prot_pilot else 'max', sec=time.time() - t0)
    if P.kind == 'tsp' and ev['feasible'].any():
        v = ev['value'][ev['feasible']]
        rec['arpd_feasible'] = float((100.0 * (v - P.meta['ref']) / P.meta['ref']).mean())
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(rec) + '\n')
    return (f"{P.name} F=min q={rec['mean_quality']:.4f} feas={rec['p_feasible']:.3f} p={rec['p_target']:.3f} | protocol F=max "
            f"q={rec['protocol_F_max']['final_mean_quality']:.4f} -> step4 {rec['step4_choice']} ({rec['sec']:.0f}s)")


def main():
    W = int(sys.argv[1]) if len(sys.argv) > 1 else 15
    jobs = [('tsp', t) for t in RB.TSP_INST] + [('gpp', g) for g in RB.GPP_INST]
    t0 = time.time()
    with mp.Pool(W) as pool:
        for k, msg in enumerate(pool.imap_unordered(job, jobs, chunksize=1)):
            if msg:
                print(f'[{k + 1}/{len(jobs)} {time.time() - t0:.0f}s] {msg}', flush=True)
    print('done', flush=True)


if __name__ == '__main__':
    main()
