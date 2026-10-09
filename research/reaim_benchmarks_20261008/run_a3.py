"""Addendum 3 runner (PROTOCOL_ADDENDUM3.md): held-out finals of the eight baselines at their papers' settings and of the
two Onsager forms at their unchanged A1 selections, fresh seeds, fresh directory results_a3/. Usage:
  python3 run_a3.py run [W]      all (problem, instance, method, budget) finals; bSB/dSB first run the papers' dt selection
                                 (64-run pilot per dt in {0.25, 0.5, 0.75, 1, 1.25}, highest pilot mean quality, ties to
                                 the smaller dt)
Seeds: finals SeedSequence([20261008, 41, p, i, m, s]); bSB/dSB dt pilots SeedSequence([20261008, 40, p, i, m, s, k]);
p = 0 Max-Cut, 1 GPP, 2 TSP; i = G-set number or 100 + TSPLIB index; m = method index in METHODS; s = budget index."""
import os

for _v in ('NUMBA_NUM_THREADS', 'OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'VECLIB_MAXIMUM_THREADS', 'MKL_NUM_THREADS'):
    os.environ[_v] = '1'

import json  # noqa: E402
import math  # noqa: E402
import multiprocessing as mp  # noqa: E402
import sys  # noqa: E402
import time  # noqa: E402
from pathlib import Path  # noqa: E402

sys.dont_write_bytecode = True
HERE = Path(__file__).resolve().parent
RES = HERE.parent
os.environ.setdefault('NUMBA_CACHE_DIR', str(HERE / '.numba_cache'))
sys.path.insert(0, str(HERE))
import numpy as np  # noqa: E402

import grids_a2 as GA  # noqa: E402  (instance lists and budgets only)
import run_a2 as RA  # noqa: E402    (problem constructors, statistics helpers)
import settings_a3 as ST  # noqa: E402
import solvers_a3 as SV3  # noqa: E402

OUT = HERE / 'results_a3'
SEED = 20261008
B_PILOT, B_FINAL = 64, 256
METHODS = GA.METHODS


def a1_selected_cfg(prob, inst, method, S):
    """The unchanged A1 selection of an Onsager form (G-set study for Max-Cut; this study's results/ otherwise)."""
    m = method.replace(' ', '_')
    if prob == 'mcp':
        p = RES / 'gset_20261007' / 'results' / 'budget' / f'G{inst}' / f'{m}_S{S}.json'
    else:
        p = HERE / 'results' / prob / RA.tag(prob, inst) / f'{m}_S{S}.json'
    return json.loads(p.read_text())['selected_cfg']


_CONST = {}


def constants(prob, inst):
    key = (prob, inst)
    if key not in _CONST:
        _CONST[key] = ST.instance_constants(RA.problem(prob, inst))
    return _CONST[key]


def stats(P, prob, ev, flips, finite, S):
    n = len(ev['feasible'])
    ok = np.ones(n, bool) if finite is None else np.asarray(finite, bool)
    feas = ev['feasible'] & ok
    q = np.where(ok, ev['quality'], 0.0)
    succ = ev['success'] & ok
    k = int(succ.sum()); p = k / n
    vals = ev['value'][feas]
    fin = dict(runs=n, n_invalid=int((~ok).sum()), ref=P.meta['ref'], target=P.meta['target'], feasible=int(feas.sum()),
               p_feasible=float(feas.mean()), p_feasible_wilson95=RA.wilson(int(feas.sum()), n), mean_quality=float(q.mean()),
               sd_quality=float(q.std(ddof=1)), mean_quality_feasible=float(q[feas].mean()) if feas.any() else None,
               k_target=k, p_target=p, p_target_wilson95=RA.wilson(k, n), k_opt=int((ev['opt'] & ok).sum()),
               best_value=(int(vals.max()) if prob == 'mcp' else int(vals.min())) if len(vals) else None, mcs99=RA.mcs99(S, p))
    if prob == 'tsp' and feas.any():
        fin['arpd_feasible'] = float((100.0 * (vals - P.meta['ref']) / P.meta['ref']).mean())
    if flips is not None:
        fin['mean_flips'] = float(np.mean(flips))
    return fin, q, feas


def job(args):
    prob, inst, method, S = args
    path = OUT / 'final' / prob / RA.tag(prob, inst) / f"{method.replace(' ', '_')}_S{S}.json"
    if path.exists():
        return None
    t0 = time.time()
    P = RA.problem(prob, inst)
    mi = METHODS.index(method); si = GA.S_LIST[prob].index(S); pc = RA.PCODE[prob]; ic = RA.icode(prob, inst)
    pilot = None
    if method in ST.OURS:
        cfg = a1_selected_cfg(prob, inst, method, S)
        source = 'A1 selection (unchanged)'
    else:
        C = constants(prob, inst)
        cands = ST.baseline_cfgs(prob, P, method, S, C)
        source = cands[0]['note']
        if len(cands) == 1:
            cfg = cands[0]
        else:
            pilot = []
            for k, c in enumerate(cands):
                s, flips, finite = SV3.run_method(P, c, B_PILOT, np.random.SeedSequence([SEED, 40, pc, ic, mi, si, k]))
                ev = P.evaluate(s)
                valid = finite is None or bool(np.all(finite))
                pilot.append(dict(dt=c['dt'], valid=valid, mean_quality=float(ev['quality'].mean()) if valid else -math.inf,
                                  p_feasible=float(ev['feasible'].mean())))
            ok = [k for k, r in enumerate(pilot) if r['valid']] or [0]
            sel = max(ok, key=lambda k: (pilot[k]['mean_quality'], -k))
            cfg = cands[sel]
    s, flips, finite = SV3.run_method(P, cfg, B_FINAL, np.random.SeedSequence([SEED, 41, pc, ic, mi, si]))
    ev = P.evaluate(s)
    fin, q, feas = stats(P, prob, ev, flips, finite, S)
    rec = dict(problem=prob, instance=P.name, method=method, S=S, cfg=cfg, source=source, dt_pilot=pilot, final=fin,
               quality=q.tolist(), values=ev['value'].tolist(), feasible_runs=feas.tolist(),
               finite=None if finite is None else np.asarray(finite).tolist(),
               flips=None if flips is None else np.asarray(flips).tolist(), sec=time.time() - t0)
    if prob == 'gpp':
        rec['imbalance'] = ev['imbalance'].tolist()
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_suffix('.tmp'); tmp.write_text(json.dumps(rec, default=list) + '\n'); tmp.rename(path)
    return (f"{prob} {RA.tag(prob, inst):6s} {method:15s} S={S:<5d} q={fin['mean_quality']:.4f} feas={fin['p_feasible']:.3f} "
            f"p={fin['p_target']:.3f} opt={fin['k_opt']} ({rec['sec']:.0f}s)")


COSTW = {'aSB': 3, 'bSB': 12, 'dSB': 12, 'ReAIM ASA': 6, 'SA': 1, 'TEC': 1.5, 'APC-SCA': 1.5}


def main():
    W = int(sys.argv[2]) if len(sys.argv) > 2 else 64
    jobs = [(prob, inst, m, S) for prob, insts in GA.INSTANCES.items() for inst in insts for m in METHODS
            for S in GA.S_LIST[prob]]
    jobs.sort(key=lambda j: -(RA.problem(j[0], j[1]).N * j[3] * COSTW.get(j[2], 1)))
    RA.run_pool(jobs, job, W, 'A3 finals')


if __name__ == '__main__':
    main()
