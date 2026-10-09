"""Addendum 5 validation V8: ReAIM ASA with ReAIM Table I's per-problem k sets (run_a5.reaim_cfg), x_best output.
  V8a  K2000 (Max-Cut, k {128, 256, 512, 1024}, T 1 -> 0.1, F = max): 4,096 and 6,400 iterations, 512 runs each
       (16 chunks of 32); P_a = P(cut >= 33,000) against ReAIM Table VII 0.47 (0.15 ms) and 0.8 (0.23 ms), the
       iteration counts as mapped by the fairness audit (time ratio 1.53, iteration ratio 1.56).
       Pass: each published P_a inside our 99 % Wilson interval or within 0.05 of ours (the audit's criterion).
  V8b  ReAIM Tables IV (MCP, 4,096), V (GPP, 4,096), VI (TSP, 8,192): 240 runs = 12 blocks of 20; best-of-20 medians
       (TSP: ARPD) against the paper. Reported; no pass/fail (the paper's numbers include ReRAM noise).
Seeds: SeedSequence([20261008, 81, crc32(job), chunk]). Results: validation_a5/<job>.json.
Usage: python3 validate_a5.py run [W] | summary"""
import os

for _v in ('NUMBA_NUM_THREADS', 'OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'VECLIB_MAXIMUM_THREADS', 'MKL_NUM_THREADS'):
    os.environ[_v] = '1'

import json  # noqa: E402
import math  # noqa: E402
import sys  # noqa: E402
import time  # noqa: E402
import zlib  # noqa: E402
from pathlib import Path  # noqa: E402

sys.dont_write_bytecode = True
HERE = Path(__file__).resolve().parent
os.environ.setdefault('NUMBA_CACHE_DIR', str(HERE / '.numba_cache'))
sys.path.insert(0, str(HERE))
import numpy as np  # noqa: E402

OUT = HERE / 'validation_a5'
T4 = ('G1', 'G2', 'G6', 'G7', 'G10', 'G11', 'G12', 'G13', 'G14', 'G19', 'G20')
T5 = (1, 2, 3, 4, 5, 14, 15, 16, 17)
T6 = ('gr17', 'gr21', 'gr24', 'fri26', 'bayg29', 'bays29')


def seed(name, chunk=0):
    return np.random.SeedSequence([20261008, 81, zlib.crc32(name.encode()), chunk])


def job(args):
    name, params = args
    path = OUT / f'{name}.json'
    if path.exists():
        return None
    t0 = time.time()
    import run_a2 as RA
    import run_a5 as R5
    import solvers_a3b as SV3B
    if params['prob'] == 'k2000':
        import k2000
        P = k2000.problem()
        cfg = R5.reaim_cfg('mcp', P, params['S'])
        rng = np.random.default_rng(seed(name, params['chunk']))
        best, final = SV3B.reaim_best(P, params['B'], cfg, rng)
        ev_b = P.evaluate(best); ev_f = P.evaluate(final)
        rec = dict(job=name, params=params, cfg=cfg, cuts_xbest=ev_b['value'].tolist(), cuts_final=ev_f['value'].tolist())
    else:
        P = RA.problem(params['prob'], params['inst'])
        cfg = R5.reaim_cfg(params['prob'], P, params['S'])
        s, _, _ = SV3B.run_method(P, cfg, params['B'], seed(name))
        ev = P.evaluate(s)
        rec = dict(job=name, params=params, cfg=cfg, values=ev['value'].tolist(), feasible=ev['feasible'].tolist())
    rec['sec'] = time.time() - t0
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_suffix('.tmp'); tmp.write_text(json.dumps(rec, default=list) + '\n'); tmp.rename(path)
    return f'{name} ({rec["sec"]:.0f}s)'


def jobs():
    J = []
    for S in (6400, 4096):
        for c in range(16):
            J.append((f'reaim_k2000_S{S}_c{c}', dict(prob='k2000', S=S, B=32, chunk=c)))
    for inst in T4:
        J.append((f'reaim_t1_t4_{inst}', dict(prob='mcp', inst=int(inst[1:]), S=4096, B=240)))
    for g in T5:
        J.append((f'reaim_t1_t5_G{g}', dict(prob='gpp', inst=g, S=4096, B=240)))
    for t in T6:
        J.append((f'reaim_t1_t6_{t}', dict(prob='tsp', inst=t, S=8192, B=240)))
    return J


def wilson(k, n, z):
    p = k / n; d = 1 + z * z / n
    c = (p + z * z / (2 * n)) / d; h = z * math.sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / d
    return c - h, c + h


def summary():
    R = {}
    for S, paper in ((4096, 0.47), (6400, 0.80)):
        xb = []; fl = []
        for c in range(16):
            d = json.loads((OUT / f'reaim_k2000_S{S}_c{c}.json').read_text()); xb += d['cuts_xbest']; fl += d['cuts_final']
        xb = np.array(xb); fl = np.array(fl); n = len(xb); k = int((xb >= 33000).sum())
        lo, hi = wilson(k, n, 2.5758)
        p = k / n
        R[f'K2000_S{S}'] = dict(n=n, p_xbest=p, wilson99=[lo, hi], p_final=float((fl >= 33000).mean()), mean_cut_xbest=float(xb.mean()),
                                paper=paper, passed=bool(lo <= paper <= hi or abs(paper - p) <= 0.05))
    print(json.dumps(R, indent=1))
    (HERE / 'validation_a5_summary.json').write_text(json.dumps(R, indent=1) + '\n')
    return R


if __name__ == '__main__':
    if sys.argv[1] == 'summary':
        summary(); sys.exit(0)
    import run_a2 as RA
    W = int(sys.argv[2]) if len(sys.argv) > 2 else 32
    RA.run_pool(jobs(), job, W, 'A5 validation V8 (ReAIM Table I k)')
