"""EXPLORATORY, POST HOC (not part of the frozen protocol): sensitivity of GPP results to the balance-penalty weight.
GPP on G1 and G14 with P in {1, 2, 8} (A = P/4; the protocol uses P = 4), S = 4096, the protocol's 32-point grids,
pilot 64 runs per point, selection by pilot mean quality, 256-run final, fresh seeds SeedSequence([20261008, 7, P, g,
method index, grid index]) and [.., 8, P, g, method index]. Methods: SA, bSB, plain SCA, TEC, APC-SCA, Onsager-kT,
Onsager-online, ReAIM. Results: explore/penalty/G<g>_P<P>/<method>.json"""
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
import problems as PB  # noqa: E402
import run_bench as RB  # noqa: E402
import solvers as SV  # noqa: E402

OUT = HERE / 'explore' / 'penalty'
METHODS = ('SA', 'bSB', 'SCA', 'TEC', 'APC-SCA', 'Onsager-kT', 'Onsager-online', 'ReAIM ASA')
S = 4096


def job(args):
    g, P_pen, method = args
    path = OUT / f'G{g}_P{P_pen}' / f"{method.replace(' ', '_')}.json"
    if path.exists():
        return None
    t0 = time.time()
    ref = json.loads((HERE / 'data' / 'gpp_reference.json').read_text())['instances'][f'G{g}']
    P = PB.gpp(g, P_pen, ref=ref['R'], target=ref['target'])
    mi = RB.METHODS.index(method)
    gr = grids.grid('gpp', method, P, S)
    pilot = []
    for gi, cfg in enumerate(gr):
        s, flips, finite = SV.run_method(P, cfg, 64, np.random.SeedSequence([20261008, 7, P_pen, g, mi, gi]))
        ev = P.evaluate(s)
        valid = finite is None or bool(np.all(finite))
        pilot.append(dict(mean_quality=float(ev['quality'].mean()) if valid else -math.inf, p_feasible=float(ev['feasible'].mean())))
    sel = max(range(len(gr)), key=lambda i: (pilot[i]['mean_quality'], -i))
    s, flips, finite = SV.run_method(P, gr[sel], 256, np.random.SeedSequence([20261008, 8, P_pen, g, mi]))
    ev = P.evaluate(s)
    rec = dict(instance=f'G{g}', P=P_pen, method=method, S=S, selected=sel, selected_cfg=gr[sel], pilot=pilot,
               mean_quality=float(ev['quality'].mean()), p_feasible=float(ev['feasible'].mean()),
               p_target=float(ev['success'].mean()), best=int(ev['value'][ev['feasible']].min()) if ev['feasible'].any() else None,
               mean_abs_imbalance=float(np.abs(ev['imbalance']).mean()), sec=time.time() - t0)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(rec) + '\n')
    return f"G{g} P={P_pen} {method:15s} q={rec['mean_quality']:.4f} feas={rec['p_feasible']:.3f} p={rec['p_target']:.3f} |M|={rec['mean_abs_imbalance']:.1f} ({rec['sec']:.0f}s)"


def main():
    W = int(sys.argv[1]) if len(sys.argv) > 1 else 32
    jobs = [(g, Pp, m) for g in (1, 14) for Pp in (1, 2, 8) for m in METHODS]
    jobs.sort(key=lambda j: -{'bSB': 4, 'ReAIM ASA': 4, 'SA': 2}.get(j[2], 1))
    t0 = time.time()
    with mp.Pool(W) as pool:
        for k, msg in enumerate(pool.imap_unordered(job, jobs, chunksize=1)):
            if msg:
                print(f'[{k + 1}/{len(jobs)} {time.time() - t0:.0f}s] {msg}', flush=True)
    print('done', flush=True)


if __name__ == '__main__':
    main()
