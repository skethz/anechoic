"""EXPLORATORY, POST HOC (after the pre-registered TEC validation failed on Fig. 3b's cycle counts): which update order or
schedule shape reproduces Du et al. 2026 Fig. 3(b) (K2000, cycles to 31,670: J_v = 0 -> 1,392; J_v = 30J -> 740)?
Variants: update order seq (index order; the pre-registered reading), rand (N single-spin updates at uniformly random
sites per cycle, as asynchronous p-bits), perm (a fresh random permutation each cycle); schedule geo or lin (100 -> 0.1
over 3,000 cycles). 64 runs per case; temporal reference = configuration at the end of the previous cycle.
Usage: python3 diag_tec.py [W]  -> diag_tec.json"""
import json
import math
import multiprocessing as mp
import os
import sys
import time
from pathlib import Path

for _v in ('NUMBA_NUM_THREADS', 'OMP_NUM_THREADS'):
    os.environ[_v] = '1'
sys.dont_write_bytecode = True
HERE = Path(__file__).resolve().parent
os.environ.setdefault('NUMBA_CACHE_DIR', str(HERE / '.numba_cache'))
sys.path.insert(0, str(HERE))
import numpy as np  # noqa: E402
import numba as nb  # noqa: E402


@nb.njit(cache=True)
def _kernel(indptr, indices, data, s, Tsched, jv, seeds, mode, W, cut_tr):
    B, N = s.shape
    S = Tsched.shape[0]
    for b in range(B):
        np.random.seed(seeds[b])
        sb = s[b].copy(); prev = sb.copy()
        h = np.zeros(N)
        for i in range(N):
            acc = 0.0
            for k in range(indptr[i], indptr[i + 1]):
                acc += data[k] * sb[indices[k]]
            h[i] = acc
        E = 0.0
        for i in range(N):
            E += -0.5 * sb[i] * h[i]
        order = np.arange(N)
        for t in range(S):
            beta = 1.0 / Tsched[t]
            if mode == 2:
                order = np.random.permutation(N)
            for n in range(N):
                if mode == 0:
                    i = n
                elif mode == 1:
                    i = np.random.randint(0, N)
                else:
                    i = order[n]
                a = 2.0 * beta * sb[i] * (h[i] + jv * prev[i])
                pf = 0.0 if a > 700.0 else (1.0 if a < -700.0 else 1.0 / (1.0 + math.exp(a)))
                if np.random.random() < pf:
                    E += 2.0 * sb[i] * h[i]
                    sb[i] = -sb[i]
                    d = 2.0 * sb[i]
                    for k in range(indptr[i], indptr[i + 1]):
                        h[indices[k]] += d * data[k]
            for i in range(N):
                prev[i] = sb[i]
            cut_tr[b, t] = (W - E) / 2.0
        s[b] = sb


def job(args):
    mode, sched, jv, chunk = args
    import k2000
    P = k2000.problem()
    rng = np.random.default_rng(np.random.SeedSequence([20261008, 70, mode, 0 if sched == 'geo' else 1, int(jv), chunk]))
    B = 8; S = 3000
    s = rng.choice(np.array([-1.0, 1.0]), size=(B, P.N))
    seeds = rng.integers(0, 2 ** 31 - 1, size=B)
    T = 100.0 * (0.001) ** (np.arange(S) / (S - 1)) if sched == 'geo' else np.linspace(100.0, 0.1, S)
    tr = np.zeros((B, S))
    W = float(P.meta['w'].sum())
    _kernel(P.indptr, P.indices, P.data.astype(np.float64), s, T, float(jv), np.asarray(seeds, np.int64), mode, W, tr)
    return (mode, sched, jv, chunk, tr.mean(0).tolist(), [int(np.argmax(c >= 31670)) + 1 if (c >= 31670).any() else None for c in tr])


def main():
    W = int(sys.argv[1]) if len(sys.argv) > 1 else 40
    jobs = [(m, sc, jv, c) for m in (0, 1, 2) for sc in ('geo', 'lin') for jv in (0.0, 30.0) for c in range(8)]
    res = {}
    t0 = time.time()
    with mp.Pool(W) as pool:
        for mode, sched, jv, chunk, tr, fp in pool.imap_unordered(job, jobs):
            key = f'{["seq", "rand", "perm"][mode]}/{sched}/jv{jv:g}'
            r = res.setdefault(key, dict(tr=[], fp=[])); r['tr'].append(tr); r['fp'] += fp
    out = {}
    for key, r in sorted(res.items()):
        m = np.mean(r['tr'], 0); fp = np.array([x if x is not None else np.inf for x in r['fp']], float)
        out[key] = dict(cross_mean_curve=int(np.argmax(m >= 31670)) + 1 if (m >= 31670).any() else None,
                        median_first_pass=float(np.median(fp)), final_mean_cut=float(m[-1]), runs=len(fp),
                        reached=int(np.isfinite(fp).sum()))
        print(key, json.dumps(out[key]), flush=True)
    (HERE / 'diag_tec.json').write_text(json.dumps(out, indent=1) + '\n')
    print(f'done {time.time() - t0:.0f}s')


if __name__ == '__main__':
    main()
