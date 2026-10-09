"""EXPLORATORY diagnostics before Addendum 4 (labelled; nothing here is a reported result):
  asb   aSB (Goto 2019) at Delta t = 0.9 and M in {2, 5}: K2000 (Nstep 186) and G-set G1, G6, G11, G14 (S = 1000):
        fraction of runs that become non-finite, the step at which |x| first exceeds 10, and the largest |x| reached;
        linear-stability numbers: lambda_max and lambda_min of J and xi0 * |lambda| * Delta t^2.
Usage: python3 diag_a4.py asb"""
import json
import math
import os
import sys
from pathlib import Path

sys.dont_write_bytecode = True
HERE = Path(__file__).resolve().parent
os.environ.setdefault('NUMBA_CACHE_DIR', str(HERE / '.numba_cache'))
sys.path.insert(0, str(HERE))
import numpy as np  # noqa: E402
import numba as nb  # noqa: E402

import k2000  # noqa: E402
import run_a2 as RA  # noqa: E402


@nb.njit(cache=True)
def _asb_diag(indptr, indices, data, x, y, S, M, delta, kick):
    B, N = x.shape
    jx = np.zeros(N, np.float32)
    one = np.float32(1.0)
    first_big = np.full(B, -1, np.int64)
    maxabs = np.zeros(B)
    for b in range(B):
        for t in range(1, S + 1):
            p = np.float32(t / S)
            for _ in range(M):
                for i in range(N):
                    x[b, i] = x[b, i] + y[b, i] * delta
                for i in range(N):
                    xv = x[b, i]
                    y[b, i] = y[b, i] - ((xv * xv + one - p) * xv) * delta
            for i in range(N):
                acc = np.float32(0.0)
                for k in range(indptr[i], indptr[i + 1]):
                    acc += data[k] * x[b, indices[k]]
                jx[i] = acc
            for i in range(N):
                y[b, i] = y[b, i] + jx[i] * kick
            m = 0.0
            for i in range(N):
                a = abs(x[b, i])
                if not np.isfinite(a):
                    a = 1e300
                if a > m:
                    m = a
            if m > maxabs[b]:
                maxabs[b] = m
            if first_big[b] < 0 and m > 10.0:
                first_big[b] = t
            if m > 1e30:
                break
    return first_big, maxabs


def asb_diag():
    out = {}
    cases = [('K2000', None, 186), ('G1', 1, 1000), ('G6', 6, 1000), ('G11', 11, 1000), ('G14', 14, 1000)]
    for name, g, S in cases:
        P = k2000.problem() if g is None else RA.problem('mcp', g)
        J = P.dense()
        ev = np.linalg.eigvalsh(J)
        sd = float(np.sqrt((J * J).mean() - J.mean() ** 2))
        xi0 = 0.7 / (sd * math.sqrt(P.N))
        for M in (2, 5):
            for dt in (0.9, 0.5):
                rng = np.random.default_rng(np.random.SeedSequence([1, 2, 3]))
                B = 16
                s0 = rng.choice(np.array([-1.0, 1.0], dtype=np.float32), size=(B, P.N))
                x = np.zeros((B, P.N), np.float32); y = (s0 * 0.1 * rng.random((B, P.N), dtype=np.float32)).astype(np.float32)
                fb, mx = _asb_diag(P.indptr, P.indices, P.data, x, y, S, M, np.float32(dt / M), np.float32(xi0 * dt))
                fin = np.isfinite(x).all(1)
                cut = P.evaluate(np.where(x >= 0, 1.0, -1.0).astype(np.float32))['value']
                rec = dict(N=P.N, S=S, M=M, dt=dt, xi0=xi0, sd_J=sd, lam_max=float(ev[-1]), lam_min=float(ev[0]),
                           kick_lmax=float(xi0 * ev[-1] * dt * dt), kick_lmin=float(xi0 * abs(ev[0]) * dt * dt),
                           frac_nonfinite=float(1 - fin.mean()), frac_exceed10=float((fb >= 0).mean()),
                           median_first_exceed10=float(np.median(fb[fb >= 0])) if (fb >= 0).any() else None,
                           median_maxabs=float(np.median(mx)), mean_cut_finite=float(cut[fin].mean()) if fin.any() else None)
                out[f'{name}/M{M}/dt{dt}'] = rec
                print(name, json.dumps(rec), flush=True)
    (HERE / 'diag_a4_asb.json').write_text(json.dumps(out, indent=1) + '\n')


if __name__ == '__main__':
    if sys.argv[1] == 'asb':
        asb_diag()
