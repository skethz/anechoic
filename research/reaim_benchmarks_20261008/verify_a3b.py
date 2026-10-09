"""Checks for Addendum 3.1 (writes verify_a3b_<host>.json):
 1. engine_sig_best runs exactly the dynamics of solvers_a3.engine_sig (identical final spins and flips) and returns the
    lowest-energy visited state: its energy equals the minimum over a numpy-reference trajectory (dense J, bias, gamma),
    including the initial state.
 2. reaim_best runs exactly the dynamics of solvers_a2.reaim_a2 (identical final state) and returns the lowest-energy
    run-phase end state, checked against an independent transcription that records every run-phase end."""
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

import problems as PB  # noqa: E402
import run_a2 as RA  # noqa: E402
import settings_a3 as ST  # noqa: E402
import solvers as SV  # noqa: E402
import solvers_a2 as SA2  # noqa: E402
import solvers_a3 as SV3  # noqa: E402
import solvers_a3b as SV3B  # noqa: E402


def ref_traj_energies(P, cfg, B, seed):
    """numpy reference of engine_sig recording E(s) for s = 1..S+1."""
    rng = np.random.default_rng(seed)
    s = rng.choice(np.array([-1.0, 1.0], dtype=np.float32), size=(B, P.N)).astype(np.float64)
    J = P.dense(); b = P.b.astype(np.float64)
    T = SV.sched(cfg['T0'], cfg['S'], cfg['tfin'])
    q = np.full(s.shape, cfg['q_reset'])
    E = lambda x: -0.5 * np.einsum('bi,ij,bj->b', x, J, x) - x @ b  # noqa: E731
    Es = [E(s)]; states = [s.copy()]
    for t, Tt in enumerate(T):
        u = rng.random(s.shape)
        z = s * (s @ J + b[None, :]) + q
        pf = 1.0 / (1.0 + np.exp(np.clip(cfg['cb'] * z / Tt, -700, 700)))
        f = u < pf
        q = np.where(f, cfg['q_reset'], np.maximum(q * cfg['r_q'], cfg['q_lim']))
        s = np.where(f, -s, s)
        Es.append(E(s)); states.append(s.copy())
    Es = np.array(Es)
    return Es.min(0), E


def main():
    out = dict(checks=[]); ok = True
    rng0 = np.random.default_rng(5)
    n = 40; Ed = [(i, j) for i in range(n) for j in range(i + 1, n) if rng0.random() < 0.2]
    gp = PB.gpp_from_edges('v_gpp40', n, np.array([e[0] for e in Ed]), np.array([e[1] for e in Ed]), np.ones(len(Ed), np.int64), 4,
                           ref=1, target=1)
    pts = rng0.random((5, 2)) * 100
    W = np.rint(np.hypot(*(pts[:, None, :] - pts[None, :, :]).transpose(2, 0, 1))).astype(np.int64); np.fill_diagonal(W, 0)
    W[(W == 0) & ~np.eye(5, dtype=bool)] = 1
    tp = PB.tsp_from_matrix('v_tsp5', W, float(W.max()), ref=1, target=1)
    mc = RA.problem('mcp', 11)
    for P, prob in ((gp, 'gpp'), (tp, 'tsp'), (mc, 'mcp')):
        C = ST.instance_constants(P)
        cfg = dict(ST.baseline_cfgs(prob, P, 'APC-SCA', 80, C)[0])
        seed = np.random.SeedSequence([9, 9])
        rng = np.random.default_rng(seed); s0 = rng.choice(np.array([-1.0, 1.0], dtype=np.float32), size=(4, P.N))
        best, f1, fin = SV3B.engine_sig_best(P, s0, cfg, rng)
        s2, f2, _ = SV3.run_method(P, cfg, 4, np.random.SeedSequence([9, 9]))
        same_dyn = bool(np.array_equal(fin, s2) and np.array_equal(f1, f2))
        emin, E = ref_traj_energies(P, cfg, 4, np.random.SeedSequence([9, 9]))
        same_best = bool(np.allclose(E(best.astype(np.float64)), emin, rtol=0, atol=1e-9))
        ok &= same_dyn and same_best
        out['checks'].append(dict(check='engine_sig_best', problem=P.name, same_dynamics=same_dyn, best_is_trajectory_min=same_best))
        print('apc best', P.name, same_dyn, same_best, flush=True)
    for prob, inst in (('gpp', 14), ('tsp', 'gr17'), ('mcp', 1)):
        P = RA.problem(prob, inst)
        C = ST.instance_constants(P)
        cfg = dict(ST.baseline_cfgs(prob, P, 'ReAIM ASA', 300, C)[0])
        xb, xf = SV3B.reaim_best(P, 3, cfg, np.random.default_rng(np.random.SeedSequence([4])))
        x2 = SA2.reaim_a2(P, 3, cfg, np.random.default_rng(np.random.SeedSequence([4])))
        same_dyn = bool(np.array_equal(xf, x2))
        # independent transcription recording run-phase ends
        rng = np.random.default_rng(np.random.SeedSequence([4]))
        Fn = {'max': np.max, 'min': np.min}[cfg['F']]
        S = cfg['S']; kset = cfg['kset']; t0, t1 = cfg['T0'], cfg['T1']; R = len(kset)
        alpha = (t1 / t0) ** (1.0 / (S - 1))
        x = rng.choice(np.array([-1.0, 1.0], dtype=np.float32), size=(3, P.N))
        cur = SV.Batch(P, x, np.full((3, 20), P.N, dtype=np.int64), Fn); t, T = 0, t0; ends = []
        while t < S:
            nn = min(32, S - t)
            rep = SV.Batch(P, np.repeat(cur.x, R, 0), np.repeat(cur.fifo, R, 0), Fn, h=np.repeat(cur.h, R, 0), M=np.repeat(cur.M, R, 0))
            kvec = np.tile(np.array(kset), 3)
            for _ in range(nn):
                rep.step(kvec, T, rng); T *= alpha
            pick = rep.energy().reshape(3, R).argmin(1); t += nn
            sel = np.arange(3) * R + pick
            cur = SV.Batch(P, rep.x[sel].copy(), rep.fifo[sel].copy(), Fn, h=rep.h[sel].copy(), M=rep.M[sel].copy())
            kb = np.array(kset)[pick]; nn = min(96, S - t)
            for _ in range(nn):
                cur.step(kb, T, rng); T *= alpha
            t += nn
            ends.append((cur.energy().copy(), cur.x.copy()))
        Es = np.array([e for e, _ in ends]); k = Es.argmin(0)
        ref_best = np.stack([ends[k[r]][1][r] for r in range(3)])
        same_best = bool(np.array_equal(xb, ref_best))
        ok &= same_dyn and same_best
        out['checks'].append(dict(check='reaim_best', problem=prob, instance=inst, same_dynamics=same_dyn, best_matches=same_best))
        print('reaim best', prob, inst, same_dyn, same_best, flush=True)
    out['all_ok'] = bool(ok)
    out['platform'] = dict(node=platform.node(), machine=platform.machine())
    (HERE / f"verify_a3b_{platform.node().split('.')[0]}.json").write_text(json.dumps(out, indent=1) + '\n')
    print('ALL OK' if ok else 'FAIL')
    return 0 if ok else 1


if __name__ == '__main__':
    sys.exit(main())
