"""Addendum 4 Table 8b re-runs (fresh seeds, fresh directory results_a4/):
  tec_seq   TEC as a sequential Glauber p-bit network (Du et al. 2026, Eq. 4), K2000 settings transferred with the
            sigma rule (OUR RULE): J_v = 30 r, T 100 r -> 0.1 r geometric over the S cycles, r = sigma/sigma_K2000.
            Seeds SeedSequence([20261008, 51, p, i, m, s]).
  asb_dt05  aSB (Goto 2019) with dt = 0.5 (the paper's stated stability bound for the explicit method; Goto 2021's aSB
            time step; OUR CHOICE), M = 2, xi0 = 0.7/(SD(J) sqrt N). Seeds [.., 52, ..].
  asb_stab  aSB with dt chosen to keep K2000's linear-stability margin (OUR RULE): dt = min(0.9, 0.9 sqrt((1 + xi0_K
            |lambda_min,K|) / (1 + xi0 |lambda_min|))), M = 2. Seeds [.., 53, ..].
Usage: python3 run_a4.py run [W]"""
import os

for _v in ('NUMBA_NUM_THREADS', 'OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'VECLIB_MAXIMUM_THREADS', 'MKL_NUM_THREADS'):
    os.environ[_v] = '1'

import json  # noqa: E402
import math  # noqa: E402
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
import solvers_a4 as SV4  # noqa: E402

OUT = HERE / 'results_a4' / 'final'
SEEDK = {'tec_seq': 51, 'asb_dt05': 52, 'asb_stab': 53}
METHOD_OF = {'tec_seq': 'TEC', 'asb_dt05': 'aSB', 'asb_stab': 'aSB'}
_K = {}


def k2000_margin():
    """1 + xi0_K |lambda_min(J_K2000)| with xi0_K = 0.7/(SD(J_K) sqrt N): K2000's stiffest-mode factor at p = 0."""
    if 'm' not in _K:
        import k2000
        C = ST.instance_constants(k2000.problem())
        xi0 = 0.7 / (C['sd_J'] * math.sqrt(C['N']))
        _K['m'] = 1.0 + xi0 * C['lam_max_minusJ']
    return _K['m']


def cfg_for(variant, prob, P, C, S):
    r = C['sigma_full'] / ST.SIGMA_K2000
    if variant == 'tec_seq':
        return dict(family='tec_seq', jv=30.0 * r, T0=100.0 * r, tfin=0.1 * r, S=S,
                    note='TEC sequential Glauber; K2000 settings x sigma/sigma_K2000 (OUR RULE)')
    xi = C['rms_J'] / C['sd_J']
    if variant == 'asb_dt05':
        return dict(family='aSB', dt=0.5, M=2, xi=xi, S=S, note='aSB dt = 0.5 (paper stability bound; OUR CHOICE), M = 2')
    xi0 = 0.7 / (C['sd_J'] * math.sqrt(P.N))
    dt = min(0.9, 0.9 * math.sqrt(k2000_margin() / (1.0 + xi0 * C['lam_max_minusJ'])))
    return dict(family='aSB', dt=dt, M=2, xi=xi, S=S, note='aSB dt keeping K2000 stability margin (OUR RULE), M = 2')


def job(args):
    variant, prob, inst, S = args
    method = METHOD_OF[variant]
    path = OUT / variant / prob / RA.tag(prob, inst) / f"{method}_S{S}.json"
    if path.exists():
        return None
    t0 = time.time()
    P = RA.problem(prob, inst)
    C = R3.constants(prob, inst)
    cfg = cfg_for(variant, prob, P, C, S)
    mi = R3.METHODS.index(method); si = GA.S_LIST[prob].index(S)
    s, flips, finite = SV4.run_method(P, cfg, R3.B_FINAL, np.random.SeedSequence([R3.SEED, SEEDK[variant], RA.PCODE[prob],
                                                                                   RA.icode(prob, inst), mi, si]))
    ev = P.evaluate(s)
    fin, q, feas = R3.stats(P, prob, ev, flips, finite, S)
    rec = dict(variant=variant, problem=prob, instance=P.name, method=method, S=S, cfg=cfg, source=cfg['note'], final=fin,
               quality=q.tolist(), values=ev['value'].tolist(), feasible_runs=feas.tolist(),
               finite=None if finite is None else np.asarray(finite).tolist(),
               flips=None if flips is None else np.asarray(flips).tolist(), sec=time.time() - t0)
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_suffix('.tmp'); tmp.write_text(json.dumps(rec, default=list) + '\n'); tmp.rename(path)
    return (f"{variant} {prob} {RA.tag(prob, inst):6s} S={S:<5d} q={fin['mean_quality']:.4f} feas={fin['p_feasible']:.3f} "
            f"p={fin['p_target']:.3f} invalid={fin['n_invalid']} ({rec['sec']:.0f}s)")


def main():
    W = int(sys.argv[2]) if len(sys.argv) > 2 else 64
    J = [(v, prob, inst, S) for v in SEEDK for prob, insts in GA.INSTANCES.items() for inst in insts for S in GA.S_LIST[prob]]
    J.sort(key=lambda j: -(RA.problem(j[1], j[2]).N * j[3] * (4 if j[0] == 'tec_seq' else 1)))
    RA.run_pool(J, job, W, 'A4 Table 8b re-runs')


if __name__ == '__main__':
    main()
