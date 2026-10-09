"""Checks for Addendum 2 before any A2 run (writes verify_a2_<host>.json):
 1. solvers_a2.reaim_a2 with F = 'max', T0 = 1.0 is bit-identical to the frozen solvers.reaim (GPP G14, TSP gr17, Max-Cut G1).
 2. The A1 points of every A2 grid, mapped by grids_a2.cfg_of, equal the frozen grids exactly and in the same order:
    grids.grid() for GPP/TSP and research/gset_20261007/run_gset.py::grid() for Max-Cut (one instance per class).
 3. The round-1 plan is deterministic (two builds give identical grids)."""
import json
import os
import platform
import sys
from pathlib import Path

sys.dont_write_bytecode = True
HERE = Path(__file__).resolve().parent
os.environ.setdefault('NUMBA_CACHE_DIR', str(HERE / '.numba_cache'))
sys.path.insert(0, str(HERE))
import numpy as np  # noqa: E402

import grids  # noqa: E402
import grids_a2 as GA  # noqa: E402
import run_a2 as RA  # noqa: E402
import solvers as SV  # noqa: E402
import solvers_a2 as SA2  # noqa: E402

NUM = ('q', 'T0', 'T1', 'tfin', 'jv', 'kappa', 'lam', 'q_reset', 'r_q', 'q_lim', 'dt', 'xi', 'S')


def same_cfg(a, b):
    if a['family'] != b['family']:
        return False
    if a['family'] == 'ReAIM':          # the frozen ReAIM had T0 = 1.0 and F = max implicitly
        a = dict(a); b = dict(b)
        for c in (a, b):
            c.setdefault('T0', 1.0); c.setdefault('F', 'max')
        if a['F'] != b['F']:
            return False
    for k in NUM:
        if (k in a) != (k in b) or (k in a and a[k] != b[k]):
            return False
    if bool(a.get('ramp', False)) != bool(b.get('ramp', False)):
        return False
    if 'kset' in a and tuple(a['kset']) != tuple(b['kset']):
        return False
    return True


def main():
    out = dict(checks=[]); ok = True
    # 1. ReAIM bit-identity
    for prob, inst, S in (('gpp', 14, 300), ('tsp', 'gr17', 300), ('mcp', 1, 300)):
        P = RA.problem(prob, inst)
        cfg = dict(family='ReAIM', kset=(1, 2, 4, 8) if prob != 'mcp' else (13, 26, 51, 102), T1=0.05, S=S)
        s1 = SV.reaim(P, 4, cfg, np.random.default_rng(np.random.SeedSequence([5, 6])))
        cfg2 = dict(cfg, F='max', T0=1.0)
        s2 = SA2.reaim_a2(P, 4, cfg2, np.random.default_rng(np.random.SeedSequence([5, 6])))
        s3, _, _ = SA2.run_method(P, cfg2, 4, np.random.SeedSequence([5, 6]))
        same = bool(np.array_equal(s1, s2) and np.array_equal(s1, s3))
        ok &= same; out['checks'].append(dict(check='reaim_a2 == solvers.reaim', problem=prob, instance=inst, identical=same))
        print('reaim', prob, inst, same, flush=True)
    # 2. A1 points == frozen grids, same order
    import run_gset as RG
    cases = [('gpp', 1, 256), ('tsp', 'gr17', 512), ('mcp', 1, 250), ('mcp', 6, 500), ('mcp', 14, 1000)]
    for prob, inst, S in cases:
        P = RA.problem(prob, inst); fam = GA.family_of(prob, inst)
        for m in GA.METHODS:
            g = GA.Grid(fam, m)
            a1 = [dict(g.fixed, **g.points[k]) for k in g.order[:g.a1_size]]
            mine = [GA.cfg_of(prob, P, m, p, S) for p in a1]
            if prob == 'mcp':
                ref = RG.grid(m, RG.graph(inst), inst, S)
            else:
                ref = grids.grid(prob, m, P, S)
            same = len(mine) == len(ref) and all(same_cfg(a, b) for a, b in zip(mine, ref))
            extra_F = all(c.get('F', 'max') == 'max' and c.get('T0', 1.0) in (1.0,) or c['family'] != 'ReAIM' for c in mine)
            ok &= bool(same and extra_F)
            out['checks'].append(dict(check='A1 points == frozen grid (values and order)', problem=prob, instance=inst,
                                      method=m, n=len(ref), identical=bool(same and extra_F)))
            print('A1', prob, inst, m, len(ref), same and extra_F, flush=True)
    # 3. deterministic plan
    g1, r1 = RA.build_grids(1); g2, r2 = RA.build_grids(1)
    same = all(json.dumps(g1[k].to_json(), default=list) == json.dumps(g2[k].to_json(), default=list) for k in g1)
    ok &= same; out['checks'].append(dict(check='round-1 plan deterministic', identical=bool(same)))
    print('plan deterministic', same)
    out['all_ok'] = bool(ok)
    import numba
    out['platform'] = dict(node=platform.node(), machine=platform.machine(), python=platform.python_version(),
                           numpy=np.__version__, numba=numba.__version__)
    (HERE / f"verify_a2_{platform.node().split('.')[0]}.json").write_text(json.dumps(out, indent=1, default=str) + '\n')
    print('ALL OK' if ok else 'FAIL')
    return 0 if ok else 1


if __name__ == '__main__':
    sys.exit(main())
