"""Reference (best-known) balanced cuts for the GPP instances G1-G5 and G14-G17, computed BEFORE any held-out run.
No published min-bisection table for these G-set instances was found; ReAIM (ISCA 2024) Table V reports best-of-20 cut
values without stating whether they are exact bisections. Reference computation (exploratory, not a benchmark result):
256 runs of SA on the GPP Ising form (P = 4) with 65,536 sweeps (16x the largest protocol budget), T 1.0 -> 0.04 sigma,
each final state repaired to exact balance (if needed) and improved by a balanced pair-swap descent (KL-style best
swap until no swap lowers the cut). R_g = the lowest balanced cut found. Writes data/gpp_reference.json."""
import hashlib
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
import numba as nb  # noqa: E402

import problems as PB  # noqa: E402
import solvers as SV  # noqa: E402

GPP = (1, 2, 3, 4, 5, 14, 15, 16, 17)
REAIM_TABLE_V = {  # ISCA 2024 Table V, best of 20 runs, per column (SFG, MFG, SFA, MFA, ASF, AMF, ASA, Neal, Tabu)
    1: (7836, 8466, 7665, 7650, 7696, 7683, 7645, 7728, 7597), 2: (7799, 8503, 7676, 7672, 7659, 7665, 7611, 7763, 7616),
    3: (7830, 8467, 7656, 7649, 7657, 7693, 7588, 7718, 7626), 4: (7849, 8481, 7668, 7683, 7673, 7651, 7635, 7707, 7620),
    5: (7806, 8471, 7671, 7683, 7676, 7697, 7664, 7751, 7612), 14: (1318, 2065, 1144, 1136, 1143, 1142, 1115, 1251, 1150),
    15: (1373, 2059, 1131, 1149, 1148, 1121, 1126, 1284, 1188), 16: (1372, 2074, 1143, 1114, 1143, 1116, 1085, 1254, 1125),
    17: (1306, 2043, 1123, 1115, 1140, 1126, 1075, 1226, 1121)}
RUNS, SWEEPS, CHUNK = 256, 65536, 8


@nb.njit(cache=True)
def balance_and_swap(s, indptr, indices, w):
    """Make s an exact bisection (move lowest-cost vertices from the larger side), then best-improvement balanced pair
    swaps until no swap reduces the cut. Returns the final cut. Gains: g_v = (cut edges at v) - (uncut edges at v)."""
    N = s.shape[0]
    g = np.zeros(N, np.int64)
    for v in range(N):
        a = 0
        for e in range(indptr[v], indptr[v + 1]):
            a += w[e] if s[indices[e]] != s[v] else -w[e]
        g[v] = a
    M = 0
    for v in range(N):
        M += s[v]
    while M != 0:
        side = 1 if M > 0 else -1
        best = -1; bg = -10 ** 9
        for v in range(N):
            if s[v] == side and g[v] > bg:
                bg = g[v]; best = v
        s[best] = -s[best]; M -= 2 * side
        g[best] = -g[best]
        for e in range(indptr[best], indptr[best + 1]):
            u = indices[e]
            g[u] += 2 * w[e] if s[u] != s[best] else -2 * w[e]
    while True:
        bgain = 0; bu = -1; bv = -1
        for u in range(N):
            if s[u] != 1:
                continue
            for v in range(N):
                if s[v] != -1:
                    continue
                wuv = 0
                for e in range(indptr[u], indptr[u + 1]):
                    if indices[e] == v:
                        wuv = w[e]
                gain = g[u] + g[v] - 2 * wuv
                if gain > bgain:
                    bgain = gain; bu = u; bv = v
        if bu < 0:
            break
        for x in (bu, bv):
            s[x] = -s[x]
            g[x] = -g[x]
            for e in range(indptr[x], indptr[x + 1]):
                u = indices[e]
                g[u] += 2 * w[e] if s[u] != s[x] else -2 * w[e]
    cut = 0
    for v in range(N):
        for e in range(indptr[v], indptr[v + 1]):
            if indices[e] > v and s[indices[e]] != s[v]:
                cut += w[e]
    return cut


def job(args):
    gi, c = args
    P = PB.gpp(gi, 4)
    a = P.sigma_T
    cfg = dict(family='SA', T0=1.0 * a, T1=0.04 * a, S=SWEEPS)
    s, _, _ = SV.run_method(P, cfg, CHUNK, np.random.SeedSequence([20261008, 99, gi, c]))
    ev = P.evaluate(s)
    # graph-only CSR (edge weights) for the swap descent
    r = np.repeat(np.arange(P.N), np.diff(P.indptr)); wv = (P.data.astype(np.int64))
    out = []
    for k in range(CHUNK):
        sk = s[k].astype(np.int64).copy()
        cut = balance_and_swap(sk, P.indptr, P.indices, wv)
        ev2 = P.evaluate(sk[None, :].astype(np.float32))
        assert ev2['feasible'][0] and int(ev2['value'][0]) == cut
        out.append(dict(raw_feasible=bool(ev['feasible'][k]), raw_cut=int(ev['cut_any'][k]), refined_cut=int(cut),
                        sha=hashlib.sha256(sk.astype(np.int8).tobytes()).hexdigest()[:16]))
    return gi, out


def main():
    W = int(sys.argv[1]) if len(sys.argv) > 1 else 64
    jobs = [(g, c) for g in GPP for c in range(RUNS // CHUNK)]
    res = {g: [] for g in GPP}; t0 = time.time()
    with mp.Pool(W) as pool:
        for gi, out in pool.imap_unordered(job, jobs, chunksize=1):
            res[gi].extend(out)
    rec = dict(method='SA (P=4, 65536 sweeps, T 1.0->0.04 sigma, 256 runs) + exact-balance repair + best balanced swap '
                      'descent', seed='SeedSequence([20261008, 99, g, chunk])', sec=time.time() - t0, instances={})
    for g in GPP:
        r = res[g]
        best = min(x['refined_cut'] for x in r)
        rec['instances'][f'G{g}'] = dict(R=best, target=int(np.floor(1.01 * best)), runs=len(r),
                                         raw_feasible_fraction=float(np.mean([x['raw_feasible'] for x in r])),
                                         raw_best_feasible=min([x['raw_cut'] for x in r if x['raw_feasible']] or [None]),
                                         refined_median=float(np.median([x['refined_cut'] for x in r])),
                                         hits_of_R=int(sum(x['refined_cut'] == best for x in r)),
                                         reaim_table_v=dict(zip(('SFG', 'MFG', 'SFA', 'MFA', 'ASF', 'AMF', 'ASA', 'Neal', 'Tabu'),
                                                                REAIM_TABLE_V[g])),
                                         reaim_table_v_min=min(REAIM_TABLE_V[g]))
        print(f"G{g}: R={best} (hits {rec['instances'][f'G{g}']['hits_of_R']}/{len(r)}), refined median "
              f"{rec['instances'][f'G{g}']['refined_median']}, raw feasible {rec['instances'][f'G{g}']['raw_feasible_fraction']:.2f}, "
              f"ReAIM Table V min {min(REAIM_TABLE_V[g])}", flush=True)
    (HERE / 'data' / 'gpp_reference.json').write_text(json.dumps(rec, indent=1) + '\n')


if __name__ == '__main__':
    main()
