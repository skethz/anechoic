"""PROTOCOL.md runner. Usage: python3 run_bench.py gpp|tsp|mcp [workers]
  gpp, tsp: all ten methods; per budget: pilot (64 runs per grid point, 32 points) -> selection by pilot mean normalized
            quality (infeasible = 0) -> 256-run held-out final. Results: results/<problem>/<instance>/<method>_S<S>.json
  mcp:      Max-Cut G1-G20 extension of research/gset_20261007 Part B: the six non-engine methods (SA, APC-SCA, ReAIM
            ASA, aSB, bSB, dSB) on the 15 instances of G1-G20 outside its 12-instance subset, with run_gset.py's own
            grid(), run_method() and stats() and its seed formulas (seed 20261007). Results: results/mcp_ext/G<g>/...
Resumable (existing result files are skipped). Single-threaded worker processes."""
import os

for _v in ('NUMBA_NUM_THREADS', 'OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'VECLIB_MAXIMUM_THREADS', 'MKL_NUM_THREADS'):
    os.environ[_v] = '1'

import itertools  # noqa: E402
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
for d in ('gset_20261007', 'ablation_20261005', 'algorithm_compare_20261007', 'reaim_reproduction_20261003',
          'theory_ideas_20261003', 'statica_reproduction_20261003'):
    if str(RES / d) not in sys.path:
        sys.path.insert(0, str(RES / d))
sys.path.insert(0, str(HERE))

import numpy as np  # noqa: E402

import grids  # noqa: E402
import problems as PB  # noqa: E402
import solvers as SV  # noqa: E402

OUT = HERE / 'results'
SEED = 20261008
B_PILOT, B_FINAL = 64, 256
METHODS = ('SA', 'SCA', 'TEC', 'APC-SCA', 'ReAIM ASA', 'aSB', 'bSB', 'dSB', 'Onsager-kT', 'Onsager-online')
ENGINE_RULES = ('SCA', 'TEC', 'Onsager-kT', 'Onsager-online')
S_LIST = {'gpp': (256, 512, 1024, 2048, 4096), 'tsp': (512, 1024, 2048, 4096, 8192)}
GPP_INST = (1, 2, 3, 4, 5, 14, 15, 16, 17)
TSP_INST = ('gr17', 'gr21', 'gr24', 'fri26', 'bayg29', 'bays29')
PCODE = {'gpp': 1, 'tsp': 2}
GPP_P = 4                                     # A = 1, B = 1 (PROTOCOL.md)


def icode(prob, inst):
    return int(inst) if prob == 'gpp' else 100 + TSP_INST.index(inst)


_PROBS = {}


def problem(prob, inst):
    key = (prob, inst)
    if key not in _PROBS:
        if prob == 'gpp':
            ref = json.loads((HERE / 'data' / 'gpp_reference.json').read_text())['instances'][f'G{inst}']
            _PROBS[key] = PB.gpp(inst, GPP_P, ref=ref['R'], target=ref['target'])
        else:
            _PROBS[key] = PB.tsp(inst)
    return _PROBS[key]


def wilson(k, n, z=1.959963984540054):
    if n == 0:
        return (0.0, 1.0)
    p = k / n; d = 1 + z * z / n
    c = (p + z * z / (2 * n)) / d; h = z * math.sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / d
    return (max(0.0, c - h), min(1.0, c + h))


def mcs99(S, p):
    return math.inf if p <= 0 else (float(S) if p >= 1 else S * math.log(0.01) / math.log1p(-p))


def stats(P, ev, flips, finite, S):
    n = len(ev['feasible'])
    ok = np.ones(n, bool) if finite is None else np.asarray(finite, bool)
    feas = ev['feasible'] & ok
    q = np.where(ok, ev['quality'], 0.0)
    succ = ev['success'] & ok
    k = int(succ.sum()); p = k / n
    vals = ev['value'][feas]
    out = dict(runs=n, n_invalid=int((~ok).sum()), ref=P.meta['ref'], target=P.meta['target'],
               feasible=int(feas.sum()), p_feasible=float(feas.mean()), p_feasible_wilson95=wilson(int(feas.sum()), n),
               mean_quality=float(q.mean()), sd_quality=float(q.std(ddof=1)),
               mean_quality_feasible=float(q[feas].mean()) if feas.any() else None,
               k_target=k, p_target=p, p_target_wilson95=wilson(k, n), k_opt=int((ev['opt'] & ok).sum()),
               best_value=(int(vals.min()) if P.kind != 'mcp' else int(vals.max())) if len(vals) else None,
               mcs99=mcs99(S, p))
    if P.kind == 'tsp' and feas.any():
        out['arpd_feasible'] = float((100.0 * (vals - P.meta['ref']) / P.meta['ref']).mean())
    if flips is not None:
        out['mean_flips'] = float(np.mean(flips))
    return out


def path_for(prob, inst, method, S):
    tag = f'G{inst}' if prob == 'gpp' else inst
    return OUT / prob / tag / f"{method.replace(' ', '_')}_S{S}.json"


def job_budget(args):
    prob, inst, method, S = args
    path = path_for(prob, inst, method, S)
    if path.exists():
        return None
    t0 = time.time()
    P = problem(prob, inst)
    mi = METHODS.index(method); si = S_LIST[prob].index(S); pc = PCODE[prob]; ic = icode(prob, inst)
    gr = grids.grid(prob, method, P, S)
    pilot = []
    for gi, cfg in enumerate(gr):
        s, flips, finite = SV.run_method(P, cfg, B_PILOT, np.random.SeedSequence([SEED, 1, pc, ic, mi, si, gi]))
        ev = P.evaluate(s)
        valid = finite is None or bool(np.all(finite))
        pilot.append(dict(cfg=cfg, valid=valid, mean_quality=float(ev['quality'].mean()) if valid else -math.inf,
                          p_feasible=float(ev['feasible'].mean()), p_target=float(ev['success'].mean()) if valid else 0.0,
                          mean_flips=float(np.mean(flips)) if flips is not None else None))
    ok = [i for i, r in enumerate(pilot) if r['valid']]
    rec = dict(problem=prob, instance=P.name, method=method, method_index=mi, S=S, N=P.N, stats=P.stats(),
               sigma_T=P.sigma_T, pilot=pilot)
    if ok:
        sel = max(ok, key=lambda i: (pilot[i]['mean_quality'], -i))
        s, flips, finite = SV.run_method(P, gr[sel], B_FINAL, np.random.SeedSequence([SEED, 2, pc, ic, mi, si]))
        ev = P.evaluate(s)
        rec.update(selected=sel, selected_cfg=gr[sel], final=stats(P, ev, flips, finite, S),
                   values=ev['value'].tolist(), feasible_runs=ev['feasible'].tolist(),
                   flips=None if flips is None else np.asarray(flips).tolist(),
                   finite=None if finite is None else np.asarray(finite).tolist())
        if prob == 'gpp':
            rec['imbalance'] = ev['imbalance'].tolist(); rec['cut_any'] = ev['cut_any'].tolist()
    else:
        rec.update(selected=None, selected_cfg=None, final=None)
    rec['sec'] = time.time() - t0
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_suffix('.tmp'); tmp.write_text(json.dumps(rec) + '\n'); tmp.rename(path)
    f = rec['final']
    if not f:
        return f'{P.name} {method} S={S}: no valid configuration ({rec["sec"]:.0f}s)'
    return (f"{P.name:7s} {method:15s} S={S:<5d} sel={sel:2d} pilotq={pilot[sel]['mean_quality']:.4f} | final q={f['mean_quality']:.4f} "
            f"feas={f['p_feasible']:.3f} p={f['p_target']:.3f} opt={f['k_opt']} best={f['best_value']} mcs99={f['mcs99']:.0f} "
            f"({rec['sec']:.0f}s)")


# ------------------------------------------------------------------ Max-Cut extension (G-set Part B procedure)
MCP_EXT = (2, 3, 4, 5, 7, 8, 9, 10, 12, 13, 15, 16, 17, 19, 20)


def job_mcp(args):
    import run_gset as RG
    g, method, S = args
    path = OUT / 'mcp_ext' / f'G{g}' / f"{method.replace(' ', '_')}_S{S}.json"
    if path.exists():
        return None
    t0 = time.time()
    G = RG.graph(g); mi = RG.METHODS.index(method); si = RG.S_LIST.index(S)
    gr = RG.grid(method, G, g, S)
    tgt = RG.BKV[f'G{g}']['target']
    pilot = []
    for gi, cfg in enumerate(gr):
        cuts, flips, finite = RG.run_method(G, cfg, RG.B_PILOT, np.random.SeedSequence([RG.SEED, 1, g, mi, si, gi]))
        valid = finite is None or bool(np.all(finite))
        pilot.append(dict(cfg=cfg, valid=valid, mean_cut=float(cuts.mean()) if valid else -math.inf,
                          p_target=float((cuts >= tgt).mean()) if valid else 0.0,
                          mean_flips=float(np.mean(flips)) if flips is not None else None))
    ok = [i for i, r in enumerate(pilot) if r['valid']]
    rec = dict(instance=f'G{g}', method=method, method_index=mi, S=S, wclass=RG.wclass(g), sigma=G.sigma, N=G.N, m=G.m,
               lam_factor=G.lam_factor, pilot=pilot, procedure='research/gset_20261007/run_gset.py job_budget (Part B)')
    if ok:
        sel = max(ok, key=lambda i: (pilot[i]['mean_cut'], -i))
        cuts, flips, finite = RG.run_method(G, gr[sel], RG.B_FINAL, np.random.SeedSequence([RG.SEED, 2, g, mi, si]))
        rec.update(selected=sel, selected_cfg=gr[sel], final=RG.stats(cuts, flips, finite, g, S), cuts=cuts.tolist(),
                   flips=None if flips is None else np.asarray(flips).tolist(),
                   finite=None if finite is None else np.asarray(finite).tolist())
    else:
        rec.update(selected=None, selected_cfg=None, final=None)
    rec['sec'] = time.time() - t0
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_suffix('.tmp'); tmp.write_text(json.dumps(rec) + '\n'); tmp.rename(path)
    f = rec['final']
    return (f"G{g:<3d} {method:15s} S={S:<5d} sel={rec['selected']} | final mean={f['mean_cut']:.1f} max={f['max_cut']} "
            f"p={f['p_target']:.3f} bkv={f['k_bkv']} mcs99={f['mcs99']:.0f} ({rec['sec']:.0f}s)" if f else
            f"G{g} {method} S={S}: no valid configuration")


def cost(prob, inst, method, S):
    n = {'gpp': 800}.get(prob) or (problem(prob, inst).N)
    w = {'aSB': 6, 'bSB': 4, 'dSB': 4, 'ReAIM ASA': 4, 'SA': 2}.get(method, 1)
    return n * S * w


def main():
    part = sys.argv[1]; W = int(sys.argv[2]) if len(sys.argv) > 2 else 32
    if part in ('gpp', 'tsp'):
        insts = GPP_INST if part == 'gpp' else TSP_INST
        jobs = [(part, i, m, S) for i in insts for m in METHODS for S in S_LIST[part]]
        jobs.sort(key=lambda j: -cost(*j))
        fn = job_budget
    elif part == 'mcp':
        import run_gset as RG
        jobs = [(g, m, S) for g in MCP_EXT for m in RG.OTHERS for S in RG.S_LIST]
        jobs.sort(key=lambda j: -(j[2] * {'aSB': 6, 'bSB': 5, 'dSB': 5, 'ReAIM ASA': 4}.get(j[1], 1)))
        fn = job_mcp
    else:
        raise SystemExit(__doc__)
    print(f'part {part}: {len(jobs)} jobs, {W} workers, start {time.strftime("%Y-%m-%d %H:%M:%S")}', flush=True)
    t0 = time.time(); done = 0
    with mp.Pool(W) as pool:
        for msg in pool.imap_unordered(fn, jobs, chunksize=1):
            done += 1
            if msg:
                print(f'[{done}/{len(jobs)} {time.time() - t0:.0f}s] {msg}', flush=True)
    print(f'part {part} complete {time.strftime("%Y-%m-%d %H:%M:%S")} ({time.time() - t0:.0f}s)', flush=True)


if __name__ == '__main__':
    main()
