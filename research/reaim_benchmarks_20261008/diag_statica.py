"""EXPLORATORY, POST HOC (after V2 gave P_a = 0.867 at the STATICA long point, against 0.827 [0.815, 0.839] in
research/statica_reproduction_20261003 with the same algorithm and graph): is the gap statistical or systematic?
  ours : solvers.engine_run family 'plain' (the engine used for STATICA/APC/Onsager), 16 chunks x 256 runs, fresh seeds
  ref  : a line-by-line transcription of statica_reproduction_20261003/statica_repro.py::sca (numpy, synchronous,
         clipped Eq. 7), 8 chunks x 256 runs, fresh seeds
Long point: K2000, S = 1,560, q = 4, T 40 -> 5 geometric, random start, final state; P_a = P(cut >= 33,000).
Usage: python3 diag_statica.py [W] -> diag_statica.json"""
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


def ref_sca(J, S, t_init, t_fin, q, sigma0, rng):
    """Transcription of statica_repro.sca (uniforms = rng.random((B, N), float32) per step)."""
    sigma = sigma0.astype(np.float32).copy()
    h = sigma @ J
    r_T = (t_fin / t_init) ** (1.0 / (S - 1))
    T = t_init
    for s in range(S):
        stay = np.clip((h * sigma + q) / (4.0 * T) + 0.5, 0.0, 1.0)
        flip = stay < rng.random(sigma.shape, dtype=np.float32)
        d = np.where(flip, -2.0 * sigma, 0.0).astype(np.float32)
        sigma += d
        h += d @ J
        T *= r_T
    return sigma


def job(args):
    which, chunk = args
    import k2000
    import solvers as SV
    P = k2000.problem()
    B = 256
    rng = np.random.default_rng(np.random.SeedSequence([20261008, 71, 0 if which == 'ours' else 1, chunk]))
    if which == 'ours':
        s0 = rng.choice(np.array([-1.0, 1.0], dtype=np.float32), size=(B, P.N))
        s, _ = SV.engine_run(P, s0, dict(family='plain', q=4.0, T0=40.0, S=1560), rng, 5.0)
    else:
        J = P.dense().astype(np.float32)
        s0 = rng.choice(np.array([-1.0, 1.0], dtype=np.float32), size=(B, P.N))
        s = ref_sca(J, 1560, 40.0, 5.0, 4.0, s0, rng)
    w = P.meta['w']; i, j = P.meta['ei'], P.meta['ej']
    cut = ((1 - s[:, i] * s[:, j]) / 2 * w).sum(1)
    return which, chunk, cut.tolist()


def main():
    W = int(sys.argv[1]) if len(sys.argv) > 1 else 24
    jobs = [('ours', c) for c in range(16)] + [('ref', c) for c in range(8)]
    out = {'ours': [], 'ref': []}
    t0 = time.time()
    with mp.Pool(W) as pool:
        for which, chunk, cut in pool.imap_unordered(job, jobs):
            out[which] += cut
    res = {}
    for k, v in out.items():
        c = np.array(v); p = float((c >= 33000).mean()); n = len(c)
        res[k] = dict(n=n, mean_cut=float(c.mean()), sd=float(c.std(ddof=1)), p=p, se_p=float(np.sqrt(p * (1 - p) / n)))
        print(k, json.dumps(res[k]), flush=True)
    (HERE / 'diag_statica.json').write_text(json.dumps(dict(summary=res, cuts=out)) + '\n')
    print(f'done {time.time() - t0:.0f}s')


if __name__ == '__main__':
    main()
