"""Addendum 3.1 runner: APC-SCA and ReAIM ASA finals with their papers' output rules (lowest-energy visited state /
x_best), same settings as Addendum 3, fresh seeds SeedSequence([20261008, 42, p, i, m, s]), fresh directory
results_a3/final_bv/. Usage: python3 run_a3b.py run [W]"""
import os

for _v in ('NUMBA_NUM_THREADS', 'OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'VECLIB_MAXIMUM_THREADS', 'MKL_NUM_THREADS'):
    os.environ[_v] = '1'

import json  # noqa: E402
import sys  # noqa: E402
import time  # noqa: E402
from pathlib import Path  # noqa: E402

sys.dont_write_bytecode = True
HERE = Path(__file__).resolve().parent
os.environ.setdefault('NUMBA_CACHE_DIR', str(HERE / '.numba_cache'))
sys.path.insert(0, str(HERE))
import numpy as np  # noqa: E402

import grids_a2 as GA  # noqa: E402
import run_a2 as RA  # noqa: E402
import run_a3 as R3  # noqa: E402
import settings_a3 as ST  # noqa: E402
import solvers_a3b as SV3B  # noqa: E402

OUT = HERE / 'results_a3' / 'final_bv'
METHODS_BV = ('APC-SCA', 'ReAIM ASA')


def job(args):
    prob, inst, method, S = args
    path = OUT / prob / RA.tag(prob, inst) / f"{method.replace(' ', '_')}_S{S}.json"
    if path.exists():
        return None
    t0 = time.time()
    P = RA.problem(prob, inst)
    mi = R3.METHODS.index(method); si = GA.S_LIST[prob].index(S)
    cfg = ST.baseline_cfgs(prob, P, method, S, R3.constants(prob, inst))[0]
    s, flips, finite = SV3B.run_method(P, cfg, R3.B_FINAL, np.random.SeedSequence([R3.SEED, 42, RA.PCODE[prob], RA.icode(prob, inst), mi, si]))
    ev = P.evaluate(s)
    fin, q, feas = R3.stats(P, prob, ev, flips, finite, S)
    rec = dict(problem=prob, instance=P.name, method=method, S=S, cfg=cfg, source=cfg['note'],
               output_rule='lowest-energy visited state (APC Alg. 2)' if method == 'APC-SCA' else 'x_best over run-phase ends (ReAIM Alg. 3)',
               final=fin, quality=q.tolist(), values=ev['value'].tolist(), feasible_runs=feas.tolist(), finite=None,
               flips=None if flips is None else np.asarray(flips).tolist(), sec=time.time() - t0)
    if prob == 'gpp':
        rec['imbalance'] = ev['imbalance'].tolist()
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_suffix('.tmp'); tmp.write_text(json.dumps(rec, default=list) + '\n'); tmp.rename(path)
    return (f"bv {prob} {RA.tag(prob, inst):6s} {method:15s} S={S:<5d} q={fin['mean_quality']:.4f} feas={fin['p_feasible']:.3f} "
            f"p={fin['p_target']:.3f} opt={fin['k_opt']} ({rec['sec']:.0f}s)")


def main():
    W = int(sys.argv[2]) if len(sys.argv) > 2 else 64
    jobs = [(prob, inst, m, S) for prob, insts in GA.INSTANCES.items() for inst in insts for m in METHODS_BV for S in GA.S_LIST[prob]]
    jobs.sort(key=lambda j: -(RA.problem(j[0], j[1]).N * j[3] * (6 if j[2] == 'ReAIM ASA' else 1.5)))
    RA.run_pool(jobs, job, W, 'A3.1 best-visited finals')


if __name__ == '__main__':
    main()
