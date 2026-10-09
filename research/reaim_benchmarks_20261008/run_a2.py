"""Addendum 2 runner (PROTOCOL_ADDENDUM2.md). Usage:
  python3 run_a2.py plan  R            build and save the round-R grids (R = 1: triggered by the A1 selections; R > 1: by
                                        the A2 union selections after round R-1); prints sizes; no runs
  python3 run_a2.py pilot R [W]        build the round-R grids (as plan) and run the 64-run pilots of every point that is new
                                        in round R (all points in round 1), fresh seeds; resumable
  python3 run_a2.py finals [W]         256-run held-out finals of the union selections after round 3, fresh seeds
  python3 run_a2.py status R           selections and edge counts after round R
Results (fresh directory): results_a2/grids_r<R>.json, results_a2/pilot/r<R>/<prob>/<inst>/<method>_S<S>_c<k>.json,
results_a2/final/<prob>/<inst>/<method>_S<S>.json. Seeds: pilots SeedSequence([20261008, 20 + R, p, i, m, s, point]),
finals SeedSequence([20261008, 31, p, i, m, s]); p = 0 Max-Cut, 1 GPP, 2 TSP; i = G-set number or 100 + TSPLIB index."""
import os

for _v in ('NUMBA_NUM_THREADS', 'OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'VECLIB_MAXIMUM_THREADS', 'MKL_NUM_THREADS'):
    os.environ[_v] = '1'

import json  # noqa: E402
import math  # noqa: E402
import multiprocessing as mp  # noqa: E402
import sys  # noqa: E402
import time  # noqa: E402
from collections import defaultdict  # noqa: E402
from pathlib import Path  # noqa: E402

sys.dont_write_bytecode = True
HERE = Path(__file__).resolve().parent
RES = HERE.parent
os.environ.setdefault('NUMBA_CACHE_DIR', str(HERE / '.numba_cache'))
sys.path.insert(0, str(HERE))

import numpy as np  # noqa: E402

import grids_a2 as GA  # noqa: E402
import problems as PB  # noqa: E402
import solvers_a2 as SA2  # noqa: E402

OUT = HERE / 'results_a2'
SEED = 20261008
B_PILOT, B_FINAL = 64, 256
PCODE = {'mcp': 0, 'gpp': 1, 'tsp': 2}
CHUNK = 8
KSET_N = {'gpp': 800, 'tsp': 841}


def icode(prob, inst):
    return 100 + GA.TSP_INST.index(inst) if prob == 'tsp' else int(inst)


def tag(prob, inst):
    return inst if prob == 'tsp' else f'G{inst}'


_PROBS = {}


def problem(prob, inst):
    key = (prob, inst)
    if key not in _PROBS:
        if prob == 'mcp':
            bkv = json.loads((RES / 'gset_20261007' / 'bkv.json').read_text())['instances']
            _PROBS[key] = PB.mcp(inst, bkv[f'G{inst}'])
        elif prob == 'gpp':
            ref = json.loads((HERE / 'data' / 'gpp_reference.json').read_text())['instances'][f'G{inst}']
            _PROBS[key] = PB.gpp(inst, 4, ref=ref['R'], target=ref['target'])
        else:
            _PROBS[key] = PB.tsp(inst)
    return _PROBS[key]


# ------------------------------------------------------------------------------------------------ A1 selections
def a1_selection(prob, inst, method, S):
    """(point, tied) of the frozen A1 selection (this study's results/, or research/gset_20261007 for Max-Cut)."""
    m = method.replace(' ', '_')
    if prob == 'mcp':
        p = RES / 'gset_20261007' / 'results' / 'budget' / f'G{inst}' / f'{m}_S{S}.json'
        if not p.exists():
            p = HERE / 'results' / 'mcp_ext' / f'G{inst}' / f'{m}_S{S}.json'
        d = json.loads(p.read_text()); key = 'mean_cut'
    else:
        d = json.loads((HERE / 'results' / prob / tag(prob, inst) / f'{m}_S{S}.json').read_text()); key = 'mean_quality'
    vals = [r[key] for r in d['pilot'] if r['valid']]
    tied = len(vals) > 1 and max(vals) == min(vals)
    return GA.point_from_cfg(prob, method, d['selected_cfg']), tied


# ------------------------------------------------------------------------------------------------ A2 selections
def pilot_records(prob, inst, method, S, upto):
    recs = {}
    for R in range(1, upto + 1):
        d = OUT / 'pilot' / f'r{R}' / prob / tag(prob, inst)
        for f in sorted(d.glob(f"{method.replace(' ', '_')}_S{S}_c*.json")):
            for r in json.loads(f.read_text())['points']:
                recs[r['idx']] = r
    return recs


def a2_selection(prob, inst, method, S, upto, grid=None):
    recs = pilot_records(prob, inst, method, S, upto)
    ok = [i for i, r in recs.items() if r['valid']]
    if grid is not None:
        assert set(recs) == set(range(grid.size())), (prob, inst, method, S, len(recs), grid.size())
    sel = max(ok, key=lambda i: (recs[i]['mean_quality'], -i))
    vals = [recs[i]['mean_quality'] for i in ok]
    tied = len(vals) > 1 and max(vals) == min(vals)
    return sel, recs[sel], tied


# ------------------------------------------------------------------------------------------------ grids per round
def build_grids(R):
    """Grids after the round-R extension (deterministic from the A1 results and the A2 pilots of rounds < R)."""
    grids = {(fam, m): GA.Grid(fam, m) for fam in GA.families() for m in GA.METHODS}
    report = {}
    for r in range(1, R + 1):
        for (fam, m), g in grids.items():
            sels = []
            for prob, inst in GA.family_instances(fam):
                for S in GA.S_LIST[prob]:
                    if r == 1:
                        sels.append(a1_selection(prob, inst, m, S))
                    else:
                        idx, rec, tied = a2_selection(prob, inst, m, S, r - 1, g)
                        sels.append((g.points[g.order[idx]], tied))
            counts = GA.edge_counts(g, [(p, t) for p, t in sels if p is not None])
            trig = GA.triggers_from(counts)
            N = 2000 if fam.startswith('mcp') else KSET_N[fam]
            applied = g.extend(trig, r, N)
            report[(r, fam, m)] = dict(counts=counts, triggers=trig, applied=applied, size=g.size())
    return grids, report


def save_grids(R, grids, report):
    OUT.mkdir(exist_ok=True)
    rec = dict(round=R, grids={f'{fam}|{m}': g.to_json() for (fam, m), g in grids.items()},
               report={f'r{r}|{fam}|{m}': v for (r, fam, m), v in report.items()})
    (OUT / f'grids_r{R}.json').write_text(json.dumps(rec, default=list) + '\n')


# ------------------------------------------------------------------------------------------------ jobs
def job_pilot(args):
    R, prob, inst, method, S, k, idxs, pts = args
    path = OUT / 'pilot' / f'r{R}' / prob / tag(prob, inst) / f"{method.replace(' ', '_')}_S{S}_c{k}.json"
    if path.exists():
        return None
    t0 = time.time()
    P = problem(prob, inst)
    mi = GA.METHODS.index(method); si = GA.S_LIST[prob].index(S)
    out = []
    for idx, p in zip(idxs, pts):
        cfg = GA.cfg_of(prob, P, method, p, S)
        s, flips, finite = SA2.run_method(P, cfg, B_PILOT, np.random.SeedSequence([SEED, 20 + R, PCODE[prob], icode(prob, inst), mi, si, idx]))
        ev = P.evaluate(s)
        valid = finite is None or bool(np.all(finite))
        out.append(dict(idx=idx, point=p, cfg=cfg, valid=valid, mean_quality=float(ev['quality'].mean()) if valid else -math.inf,
                        p_feasible=float(ev['feasible'].mean()), p_target=float(ev['success'].mean()) if valid else 0.0,
                        mean_flips=float(np.mean(flips)) if flips is not None else None))
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_suffix('.tmp'); tmp.write_text(json.dumps(dict(round=R, points=out, sec=time.time() - t0), default=list) + '\n')
    tmp.rename(path)
    return f'r{R} {prob} {tag(prob, inst):6s} {method:15s} S={S:<5d} chunk {k}: {len(idxs)} pts ({time.time() - t0:.0f}s)'


def wilson(k, n, z=1.959963984540054):
    if n == 0:
        return (0.0, 1.0)
    p = k / n; d = 1 + z * z / n
    c = (p + z * z / (2 * n)) / d; h = z * math.sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / d
    return (max(0.0, c - h), min(1.0, c + h))


def mcs99(S, p):
    return math.inf if p <= 0 else (float(S) if p >= 1 else S * math.log(0.01) / math.log1p(-p))


def job_final(args):
    prob, inst, method, S, sel, cfg, pilot_q = args
    path = OUT / 'final' / prob / tag(prob, inst) / f"{method.replace(' ', '_')}_S{S}.json"
    if path.exists():
        return None
    t0 = time.time()
    P = problem(prob, inst)
    mi = GA.METHODS.index(method); si = GA.S_LIST[prob].index(S)
    s, flips, finite = SA2.run_method(P, cfg, B_FINAL, np.random.SeedSequence([SEED, 31, PCODE[prob], icode(prob, inst), mi, si]))
    ev = P.evaluate(s)
    n = len(ev['feasible'])
    ok = np.ones(n, bool) if finite is None else np.asarray(finite, bool)
    feas = ev['feasible'] & ok
    q = np.where(ok, ev['quality'], 0.0)
    succ = ev['success'] & ok
    k = int(succ.sum()); p = k / n
    vals = ev['value'][feas]
    fin = dict(runs=n, n_invalid=int((~ok).sum()), ref=P.meta['ref'], target=P.meta['target'], feasible=int(feas.sum()),
               p_feasible=float(feas.mean()), p_feasible_wilson95=wilson(int(feas.sum()), n), mean_quality=float(q.mean()),
               sd_quality=float(q.std(ddof=1)), mean_quality_feasible=float(q[feas].mean()) if feas.any() else None,
               k_target=k, p_target=p, p_target_wilson95=wilson(k, n), k_opt=int((ev['opt'] & ok).sum()),
               best_value=(int(vals.max()) if prob == 'mcp' else int(vals.min())) if len(vals) else None, mcs99=mcs99(S, p))
    if prob == 'tsp' and feas.any():
        fin['arpd_feasible'] = float((100.0 * (vals - P.meta['ref']) / P.meta['ref']).mean())
    if flips is not None:
        fin['mean_flips'] = float(np.mean(flips))
    rec = dict(problem=prob, instance=P.name, method=method, S=S, selected=sel, selected_cfg=cfg, pilot_mean_quality=pilot_q,
               final=fin, quality=q.tolist(), values=ev['value'].tolist(), feasible_runs=feas.tolist(),
               flips=None if flips is None else np.asarray(flips).tolist(), sec=time.time() - t0)
    if prob == 'gpp':
        rec['imbalance'] = ev['imbalance'].tolist()
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_suffix('.tmp'); tmp.write_text(json.dumps(rec, default=list) + '\n'); tmp.rename(path)
    return (f"final {prob} {tag(prob, inst):6s} {method:15s} S={S:<5d} sel={sel} q={fin['mean_quality']:.4f} "
            f"feas={fin['p_feasible']:.3f} p={fin['p_target']:.3f} ({rec['sec']:.0f}s)")


COSTW = {'aSB': 6, 'bSB': 5, 'dSB': 5, 'ReAIM ASA': 4, 'TEC': 2, 'Onsager-kT': 2, 'SA': 1}


def cost(prob, inst, method, S, n):
    N = problem(prob, inst).N
    return N * S * n * COSTW.get(method, 1.5)


def run_pool(jobs, fn, W, label):
    print(f'{label}: {len(jobs)} jobs, {W} workers, start {time.strftime("%Y-%m-%d %H:%M:%S")}', flush=True)
    t0 = time.time(); done = 0
    with mp.Pool(W) as pool:
        for msg in pool.imap_unordered(fn, jobs, chunksize=1):
            done += 1
            if msg:
                print(f'[{done}/{len(jobs)} {time.time() - t0:.0f}s] {msg}', flush=True)
    print(f'{label} complete {time.strftime("%Y-%m-%d %H:%M:%S")} ({time.time() - t0:.0f}s)', flush=True)


def plan(R):
    grids, report = build_grids(R)
    save_grids(R, grids, report)
    tot = 0
    for (fam, m), g in grids.items():
        prev = None
        line = f'{fam:8s} {m:15s} A1={g.a1_size:3d} cap={g.cap():3d} size={g.size():3d}'
        adds = [f"r{e['round']}:{e['axis']}-{e['end']}={e.get('value')}({e['result']})" for e in g.log]
        print(line, ' '.join(adds))
        tot += g.size() * len(GA.family_instances(fam)) * 5
    print('total point-evaluations (all rounds, all instances/budgets):', tot)
    return grids


def pilot(R, W):
    grids, report = build_grids(R)
    save_grids(R, grids, report)
    prev = build_grids(R - 1)[0] if R > 1 else None
    jobs = []
    for (fam, m), g in grids.items():
        start = 0 if R == 1 else prev[(fam, m)].size()
        new = list(range(start, g.size()))
        if not new:
            continue
        for prob, inst in GA.family_instances(fam):
            for S in GA.S_LIST[prob]:
                for k in range(0, len(new), CHUNK):
                    ii = new[k:k + CHUNK]
                    pts = [dict(g.fixed, **g.points[g.order[i]]) for i in ii]
                    jobs.append((R, prob, inst, m, S, k // CHUNK, ii, pts))
    jobs.sort(key=lambda j: -cost(j[1], j[2], j[3], j[4], len(j[6])))
    run_pool(jobs, job_pilot, W, f'A2 pilot round {R}')


def finals(W):
    grids, _ = build_grids(GA.MAX_ROUNDS)
    jobs = []
    for (fam, m), g in grids.items():
        for prob, inst in GA.family_instances(fam):
            P = problem(prob, inst)
            for S in GA.S_LIST[prob]:
                sel, rec, tied = a2_selection(prob, inst, m, S, GA.MAX_ROUNDS, g)
                jobs.append((prob, inst, m, S, sel, rec['cfg'], rec['mean_quality']))
    jobs.sort(key=lambda j: -cost(j[0], j[1], j[2], j[3], 4))
    run_pool(jobs, job_final, W, 'A2 finals')


def status(R):
    grids, report = build_grids(R)
    for (fam, m), g in grids.items():
        sels = []
        for prob, inst in GA.family_instances(fam):
            for S in GA.S_LIST[prob]:
                idx, rec, tied = a2_selection(prob, inst, m, S, R, g)
                sels.append((g.points[g.order[idx]], tied))
        c = GA.edge_counts(g, sels)
        print(f'{fam:8s} {m:15s} size={g.size():3d} ' + '; '.join(f"{k}: {v['min']}/{v['n']} min, {v['max']}/{v['n']} max" for k, v in c.items()))


if __name__ == '__main__':
    cmd = sys.argv[1]
    if cmd == 'plan':
        plan(int(sys.argv[2]))
    elif cmd == 'pilot':
        pilot(int(sys.argv[2]), int(sys.argv[3]) if len(sys.argv) > 3 else 64)
    elif cmd == 'finals':
        finals(int(sys.argv[2]) if len(sys.argv) > 2 else 64)
    elif cmd == 'status':
        status(int(sys.argv[2]))
    else:
        raise SystemExit(__doc__)
