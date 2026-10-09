"""EXPLORATORY diagnostic (not pre-registered, not used for any number in a table): which reading of TEC's unstated
details reproduces its Fig. 3(b) (cycles for the mean cut to reach 31,670: J_v = 0: ~1,392; J_v = 30: ~740)?
Variants of the sequential Glauber kernel, K2000, 3,000 cycles, k_BT 100 -> 0.1, 32 runs per (variant, J_v):
  geom     geometric T (PROTOCOL_PUB reading)
  lin      linear T
  geom_b1  geometric T, flip probability 1/(1 + exp{beta s (h + J_v s_prev)}) (no factor 2 in the exponent)
  lin_b1   linear T, no factor 2
  geom_a05 geometric T, each spin attempts an update with probability 0.5 per cycle (switching rate alpha)
Seeds SeedSequence([999, 7, variant, J_v index]) (outside every pre-registered seed space)."""
import json
import math
import sys
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import numpy as np
import numba as nb

sys.path.insert(0, str(Path(__file__).resolve().parent))
import pub_methods as P  # noqa: E402
from common import M  # noqa: E402

HERE = Path(__file__).resolve().parent
VARIANTS = ('geom', 'lin', 'geom_b1', 'lin_b1', 'geom_a05')
JVS = (0.0, 30.0)
S, B = 3000, 32


@nb.njit(cache=True)
def kern(J, s, Tsched, Jv, seeds, E, fac, attempt):
    B, N = s.shape
    for b in range(B):
        np.random.seed(seeds[b])
        sb = s[b].copy(); h = np.zeros(N)
        for i in range(N):
            acc = 0.0
            for j in range(N):
                acc += J[i, j] * sb[j]
            h[i] = acc
        prev = sb.copy()
        for t in range(Tsched.shape[0]):
            beta = 1.0 / Tsched[t]
            for i in range(N):
                prev[i] = sb[i]
            for i in range(N):
                if attempt < 1.0 and np.random.random() >= attempt:
                    continue
                x = fac * beta * sb[i] * (h[i] + Jv * prev[i])
                pf = 0.0 if x > 700.0 else (1.0 if x < -700.0 else 1.0 / (1.0 + math.exp(x)))
                if np.random.random() < pf:
                    sb[i] = -sb[i]; d = 2.0 * sb[i]
                    for j in range(N):
                        h[j] += d * J[i, j]
            e = 0.0
            for i in range(N):
                e += sb[i] * h[i]
            E[t + 1, b] = -0.5 * e


def job(arg):
    vi, ji = arg
    v = VARIANTS[vi]; J, sumw = M.m.load()
    s0, seeds = P.draw_sequential_inputs(np.random.SeedSequence([999, 7, vi, ji]), B, len(J))
    Ts = P.geom_sched(100.0, 0.1, S) if v.startswith('geom') else np.linspace(100.0, 0.1, S)
    fac = 1.0 if v.endswith('b1') else 2.0
    att = 0.5 if v.endswith('a05') else 1.0
    s = s0.astype(np.float64).copy(); E = np.zeros((S + 1, B))
    kern(J.astype(np.float64), s, Ts, JVS[ji], np.asarray(seeds, np.int64), E, fac, att)
    c = (sumw - E[1:]) / 2
    mean = c.mean(1)
    first = int(np.argmax(mean >= 31670)) + 1 if (mean >= 31670).any() else None
    return v, JVS[ji], first, float(mean[-1])


def main(workers):
    out = {}
    with ProcessPoolExecutor(max_workers=workers) as ex:
        for v, jv, first, fin in ex.map(job, [(vi, ji) for vi in range(len(VARIANTS)) for ji in range(len(JVS))]):
            out.setdefault(v, {})[str(jv)] = dict(first_cycle=first, final_mean_cut=fin)
            print(v, jv, first, fin, flush=True)
    for v in out:
        a, b = out[v]['0.0']['first_cycle'], out[v]['30.0']['first_cycle']
        out[v]['ratio'] = a / b if a and b else None
    (HERE / 'results_validation').mkdir(exist_ok=True)
    (HERE / 'results_validation' / 'diag_tec_schedule_EXPLORATORY.json').write_text(json.dumps(out, indent=1) + '\n')
    print(json.dumps(out))


if __name__ == '__main__':
    main(int(sys.argv[1]) if len(sys.argv) > 1 else 10)
