"""Checks for Addendum 3 before any A3 run (writes verify_a3_<host>.json):
 1. solvers_a3.engine_sig (APC-SCA and TEC, exact logistic) equals a plain numpy reference with a dense J, bias and the
    uniform GPP coupling, bit for bit in spins and flips (same RNG draws), on small GPP, TSP and Max-Cut cases.
 2. settings_a3.instance_constants' Neal beta range equals an independent dictionary-based transcription of
    dwave-samplers 1.2.0 `_default_ising_beta_range` (Ising h = -b, J_ij = -J_ij for i < j) on GPP G14, TSP gr17, G1.
 3. reaim_a2 with F = min, T0 = 0.5 equals an independent transcription (explore_reaim_F.reaim_F with T0 inserted).
 4. One short run of every baseline and Onsager form on one instance per problem (no exceptions, finite results)."""
import json
import math
import os
import platform
import sys
from collections import defaultdict
from pathlib import Path

sys.dont_write_bytecode = True
HERE = Path(__file__).resolve().parent
os.environ.setdefault('NUMBA_CACHE_DIR', str(HERE / '.numba_cache'))
sys.path.insert(0, str(HERE))
import numpy as np  # noqa: E402

import problems as PB  # noqa: E402
import run_a2 as RA  # noqa: E402
import run_a3 as R3  # noqa: E402
import settings_a3 as ST  # noqa: E402
import solvers as SV  # noqa: E402
import solvers_a2 as SA2  # noqa: E402
import solvers_a3 as SV3  # noqa: E402


def ref_engine_sig(P, cfg, B, seed):
    rng = np.random.default_rng(seed)
    s = rng.choice(np.array([-1.0, 1.0], dtype=np.float32), size=(B, P.N)).astype(np.float64)
    J = P.dense(); b = P.b.astype(np.float64)
    S = cfg['S']; T = SV.sched(cfg['T0'], S, cfg['tfin'])
    apc = cfg['family'] == 'apc_sig'; tec = cfg['family'] == 'tec_sig'
    q = np.full(s.shape, cfg['q_reset'] if apc else cfg.get('q', 0.0))
    sp = s.copy(); flips = np.zeros(B)
    for t, Tt in enumerate(T):
        u = rng.random(s.shape)
        fld = s @ J + b[None, :]
        if tec and t > 0:
            fld = fld + cfg['jv'] * sp
        z = s * fld + q
        a = np.clip(cfg['cb'] * z / Tt, -700, 700)
        pf = 1.0 / (1.0 + np.exp(a))
        f = u < pf
        if apc:
            q = np.where(f, cfg['q_reset'], np.maximum(q * cfg['r_q'], cfg['q_lim']))
        sp = s.copy()
        s = np.where(f, -s, s); flips += f.sum(1)
    return s.astype(np.float32), flips


def neal_ref(P):
    """dwave-samplers 1.2.0 _default_ising_beta_range on the Ising model (h, J) of P (dictionary transcription)."""
    J = P.dense(); N = P.N
    h = {i: -float(P.b[i]) for i in range(N)}
    Jd = {(i, j): -float(J[i, j]) for i in range(N) for j in range(i + 1, N) if J[i, j] != 0}
    sum_abs_bias_dict = defaultdict(int, {k: abs(v) for k, v in h.items()})
    min_abs_bias_dict = {k: v for k, v in sum_abs_bias_dict.items() if v != 0}
    for (k1, k2), v in Jd.items():
        for k in (k1, k2):
            sum_abs_bias_dict[k] += abs(v)
            if v != 0:
                min_abs_bias_dict[k] = min(abs(v), min_abs_bias_dict[k]) if k in min_abs_bias_dict else abs(v)
    max_effective_field = max(sum_abs_bias_dict.values(), default=0)
    hot_beta = np.log(2) / (2 * max_effective_field)
    values_array = np.array(list(min_abs_bias_dict.values()), dtype=float)
    min_effective_field = np.min(values_array)
    number_min_gaps = np.sum(min_effective_field == values_array)
    cold_beta = np.log(number_min_gaps / 0.01) / (2 * min_effective_field)
    return float(hot_beta), float(cold_beta)


def reaim_F_ref(P, B, cfg, rng):
    """Independent transcription: explore_reaim_F.reaim_F with the initial temperature T0 of cfg."""
    F = {'max': np.max, 'min': np.min}[cfg['F']]
    S = cfg['S']; kset = cfg['kset']; t0, t1 = cfg['T0'], cfg['T1']; it_trial, it_run = 32, 96
    N = P.N; R = len(kset)
    alpha = (t1 / t0) ** (1.0 / (S - 1))
    x = rng.choice(np.array([-1.0, 1.0], dtype=np.float32), size=(B, N))
    cur = SV.Batch(P, x, np.full((B, 20), N, dtype=np.int64), F)
    t, T = 0, t0
    while t < S:
        n = min(it_trial, S - t)
        rep = SV.Batch(P, np.repeat(cur.x, R, 0), np.repeat(cur.fifo, R, 0), F, h=np.repeat(cur.h, R, 0), M=np.repeat(cur.M, R, 0))
        kvec = np.tile(np.array(kset), B)
        for _ in range(n):
            rep.step(kvec, T, rng); T *= alpha
        pick = rep.energy().reshape(B, R).argmin(1); t += n
        sel = np.arange(B) * R + pick
        cur = SV.Batch(P, rep.x[sel].copy(), rep.fifo[sel].copy(), F, h=rep.h[sel].copy(), M=rep.M[sel].copy())
        kbest = np.array(kset)[pick]
        n = min(it_run, S - t)
        for _ in range(n):
            cur.step(kbest, T, rng); T *= alpha
        t += n
    return cur.x


def main():
    out = dict(checks=[]); ok = True
    rng0 = np.random.default_rng(3)
    n = 40; E = [(i, j) for i in range(n) for j in range(i + 1, n) if rng0.random() < 0.2]
    gp = PB.gpp_from_edges('v_gpp40', n, np.array([e[0] for e in E]), np.array([e[1] for e in E]), np.ones(len(E), np.int64), 4,
                           ref=1, target=1)
    pts = rng0.random((5, 2)) * 100
    W = np.rint(np.hypot(*(pts[:, None, :] - pts[None, :, :]).transpose(2, 0, 1))).astype(np.int64); np.fill_diagonal(W, 0)
    W[(W == 0) & ~np.eye(5, dtype=bool)] = 1
    tp = PB.tsp_from_matrix('v_tsp5', W, float(W.max()), ref=1, target=1)
    mc = RA.problem('mcp', 11)
    for P in (gp, tp, mc):
        C = ST.instance_constants(P)
        for m in ('APC-SCA', 'TEC'):
            for cfg in ST.baseline_cfgs('tsp' if P is tp else ('gpp' if P is gp else 'mcp'), P, m, 60, C):
                s1, f1, _ = SV3.run_method(P, cfg, 3, np.random.SeedSequence([7, 7]))
                s2, f2 = ref_engine_sig(P, cfg, 3, np.random.SeedSequence([7, 7]))
                same = bool(np.array_equal(s1, s2) and np.array_equal(f1, f2))
                ok &= same; out['checks'].append(dict(check='engine_sig == numpy reference', problem=P.name, method=m, identical=same,
                                                      flips=float(np.mean(f1))))
                print('engine_sig', P.name, m, same, float(np.mean(f1)), flush=True)
    for prob, inst in (('gpp', 14), ('tsp', 'gr17'), ('mcp', 1)):
        P = RA.problem(prob, inst)
        C = ST.instance_constants(P)
        hb, cb = neal_ref(P)
        same = math.isclose(hb, C['neal']['hot_beta'], rel_tol=1e-12) and math.isclose(cb, C['neal']['cold_beta'], rel_tol=1e-12)
        ok &= same; out['checks'].append(dict(check='Neal default beta range == transcription', problem=prob, instance=inst,
                                              identical=bool(same), hot=hb, cold=cb))
        print('neal', prob, inst, same, hb, cb, flush=True)
    for prob, inst in (('gpp', 14), ('tsp', 'gr17')):
        P = RA.problem(prob, inst)
        cfg = dict(family='ReAIM', kset=(1, 2, 6, 16), T0=0.5, T1=0.1, F='min', S=200)
        s1 = SA2.reaim_a2(P, 3, cfg, np.random.default_rng(np.random.SeedSequence([8])))
        s2 = reaim_F_ref(P, 3, cfg, np.random.default_rng(np.random.SeedSequence([8])))
        same = bool(np.array_equal(s1, s2)); ok &= same
        out['checks'].append(dict(check='reaim_a2(F=min, T0=0.5) == transcription', problem=prob, instance=inst, identical=same))
        print('reaim', prob, inst, same, flush=True)
    # 4. smoke: every method, one instance per problem, small B
    R3.B_PILOT, R3.B_FINAL = 4, 8
    import tempfile
    R3.OUT = Path(tempfile.mkdtemp(dir=os.environ.get('TMPDIR')))
    for prob, inst, S in (('mcp', 14, 250), ('gpp', 14, 256), ('tsp', 'gr17', 512)):
        for m in R3.METHODS:
            msg = R3.job((prob, inst, m, S))
            print('smoke', msg, flush=True)
            out['checks'].append(dict(check='smoke run', problem=prob, method=m, msg=msg))
    out['all_ok'] = bool(ok)
    import numba
    out['platform'] = dict(node=platform.node(), machine=platform.machine(), numpy=np.__version__, numba=numba.__version__)
    (HERE / f"verify_a3_{platform.node().split('.')[0]}.json").write_text(json.dumps(out, indent=1, default=str) + '\n')
    print('ALL OK' if ok else 'FAIL')
    return 0 if ok else 1


if __name__ == '__main__':
    sys.exit(main())
