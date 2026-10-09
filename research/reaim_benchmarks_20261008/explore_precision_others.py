"""EXPLORATORY, POST HOC (not part of the frozen protocol): coupling precision needed by the problems themselves, measured
with the two methods that solve GPP and TSP at full precision (SA and bSB; float model, solvers.py). For each GPP and TSP
instance, the held-out-selected configuration at ReAIM's budget (GPP S = 4096, TSP S = 8192) is run on the K-bit
quantized instance (problems.Problem.quantize, as in the precision study) for K in {2, 3, 4, 6, 8} and full precision:
SA temperatures are multiplied by alpha_K; bSB is scale-free (Goto normalisation by sigma_J of the quantized J).
256 runs per K with the same seed for every K, SeedSequence([20261008, 11, p, i, method index]). The objective is
evaluated on the true problem. Results: explore/precision_others/<problem>_<instance>_<method>.json"""
import json
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

import run_bench as RB  # noqa: E402
import solvers as SV  # noqa: E402

OUT = HERE / 'explore' / 'precision_others'
KS = (None, 2, 3, 4, 6, 8)
S_OF = {'gpp': 4096, 'tsp': 8192}


def job(args):
    prob, inst, method = args
    tag = f'G{inst}' if prob == 'gpp' else inst
    path = OUT / f'{prob}_{tag}_{method}.json'
    if path.exists():
        return None
    t0 = time.time()
    P = RB.problem(prob, inst); S = S_OF[prob]
    d = json.loads((HERE / 'results' / prob / tag / f'{method}_S{S}.json').read_text())
    cfg0 = d['selected_cfg']
    res = {}
    for K in KS:
        Q, alpha = P.quantize(K)
        cfg = dict(cfg0)
        if method == 'SA':
            cfg['T0'] = cfg0['T0'] * alpha; cfg['T1'] = cfg0['T1'] * alpha
        s, _, _ = SV.run_method(Q, cfg, 256, np.random.SeedSequence([20261008, 11, RB.PCODE[prob], RB.icode(prob, inst),
                                                                      RB.METHODS.index(method)]))
        ev = P.evaluate(s)
        res[str(K)] = dict(alpha=alpha, mean_quality=float(ev['quality'].mean()), p_feasible=float(ev['feasible'].mean()),
                           p_target=float(ev['success'].mean()), quality=ev['quality'].tolist())
    rec = dict(problem=prob, instance=P.name, method=method, S=S, selected_cfg=cfg0, by_K=res, sec=time.time() - t0)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(rec) + '\n')
    return f"{P.name} {method}: " + ' '.join(f"K={k}:{v['mean_quality']:.3f}/{v['p_feasible']:.2f}" for k, v in res.items()) + f" ({rec['sec']:.0f}s)"


def main():
    W = int(sys.argv[1]) if len(sys.argv) > 1 else 16
    jobs = [('tsp', t, m) for t in RB.TSP_INST for m in ('bSB', 'SA')] + [('gpp', g, m) for g in RB.GPP_INST for m in ('bSB', 'SA')]
    t0 = time.time()
    with mp.Pool(W) as pool:
        for k, msg in enumerate(pool.imap_unordered(job, jobs, chunksize=1)):
            if msg:
                print(f'[{k + 1}/{len(jobs)} {time.time() - t0:.0f}s] {msg}', flush=True)
    print('done', flush=True)


if __name__ == '__main__':
    main()
