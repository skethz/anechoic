"""Verification before any benchmark run (writes verify.json).
A. Formulations by brute force over all 2^N states: GPP (N <= 14) and TSP (n = 3, 4; N = 9, 16) energies equal the
   Lucas objective up to a constant; fields equal the dense J s + b; ground states are optimal feasible solutions.
B. b = 0, gamma = 0: every kernel is bit-identical to research/gset_20261007 engine.run / others.* (G-set G1, G14 and
   synthetic graphs): engine families plain, tec, tecT, onsager, apc; SA; ReAIM; dSB; bSB; aSB.
C. Bias and uniform coupling: engine_run equals a dense copy of abl.run with h(0) = s @ J + b (bit for bit); SA, ReAIM and
   dSB equal dense references (legacy MT19937 stream / dense float32 products / exact integer sums); every family is
   checked on small GPP (gamma != 0) and TSP (bias != 0) instances.
"""
import json
import math
import os
import sys
import time
from pathlib import Path

sys.dont_write_bytecode = True
HERE = Path(__file__).resolve().parent
RES = HERE.parent
os.environ.setdefault('NUMBA_CACHE_DIR', str(HERE / '.numba_cache'))
for d in ('gset_20261007', 'ablation_20261005', 'algorithm_compare_20261007', 'reaim_reproduction_20261003',
          'theory_ideas_20261003', 'statica_reproduction_20261003'):
    if str(RES / d) not in sys.path:
        sys.path.insert(0, str(RES / d))
sys.path.insert(0, str(HERE))

import itertools  # noqa: E402

import numpy as np  # noqa: E402

import problems as PB  # noqa: E402
import solvers as SV  # noqa: E402

OUT = {}


def log(*a):
    print(*a, flush=True)


# ============================================================ A. formulations
def all_states(N):
    k = np.arange(2 ** N)[:, None]
    return np.where((k >> np.arange(N)[None, :]) & 1, 1.0, -1.0)


def brute_gpp(n, p, P, seed):
    rng = np.random.default_rng(seed)
    E = [(i, j) for i in range(n) for j in range(i + 1, n) if rng.random() < p]
    i = np.array([e[0] for e in E]); j = np.array([e[1] for e in E]); w = np.ones(len(E), np.int64)
    pr = PB.gpp_from_edges(f'bf_gpp{n}', n, i, j, w, P)
    S = all_states(n)
    e = pr.energy(S)
    J = pr.dense()
    e2 = -0.5 * np.einsum('bi,ij,bj->b', S, J, S)
    cut = ((S[:, i] != S[:, j]).astype(int)).sum(1)
    M = S.sum(1)
    const = e - (2 * cut + 0.5 * P * M ** 2)
    feas = M == 0
    minbis = cut[feas].min()
    gs = np.flatnonzero(e == e.min())
    h_ok = np.allclose(pr.field(S[:64]), S[:64] @ J)
    res = dict(n=n, m=len(E), P=P, energy_matches_dense=bool(np.allclose(e, e2)),
               energy_minus_objective_constant=bool(np.ptp(const) == 0), field_ok=bool(h_ok),
               lucas_bound_P=int(max(np.bincount(np.r_[i, j], minlength=n).max(), 1)),
               min_bisection=int(minbis), ground_states_balanced=bool(np.all(M[gs] == 0)),
               ground_state_cut=int(cut[gs].min()))
    res['ok'] = bool(res['energy_matches_dense'] and res['energy_minus_objective_constant'] and res['field_ok']
                     and (P < res['lucas_bound_P'] or (res['ground_states_balanced'] and res['ground_state_cut'] == minbis)))
    return res


def brute_tsp(n, seed, A_rule):
    rng = np.random.default_rng(seed)
    pts = rng.random((n, 2)) * 100
    W = np.rint(np.hypot(*(pts[:, None, :] - pts[None, :, :]).transpose(2, 0, 1))).astype(np.int64)
    W[W == 0] = 1; np.fill_diagonal(W, 0)
    A = float(W.max()) * A_rule
    pr = PB.tsp_from_matrix(f'bf_tsp{n}', W, A, ref=1, target=1)
    N = n * n
    S = all_states(N)
    e = pr.energy(S)
    J = pr.dense(); b = pr.b.astype(np.float64)
    e2 = -0.5 * np.einsum('bi,ij,bj->b', S, J, S) - S @ b
    x = (S.reshape(-1, n, n) > 0).astype(np.int64)
    Hpen = A * (((1 - x.sum(2)) ** 2).sum(1) + ((1 - x.sum(1)) ** 2).sum(1))
    Hd = np.zeros(len(S))
    for j in range(n):
        Hd += np.einsum('bu,uv,bv->b', x[:, :, j], W, x[:, :, (j + 1) % n])
    H = Hpen + Hd
    const = e - 4 * H
    feas = np.all(x.sum(1) == 1, axis=1) & np.all(x.sum(2) == 1, axis=1)
    perms = list(itertools.permutations(range(n)))
    Lopt = min(sum(W[p[k], p[(k + 1) % n]] for k in range(n)) for p in perms)
    emin_feas = e[feas].min(); emin_all = e.min()
    gs = np.flatnonzero(e == emin_all)
    res = dict(n=n, A=A, A_over_maxW=A_rule, energy_matches_dense=bool(np.allclose(e, e2)),
               energy_minus_4H_constant=bool(np.ptp(const) == 0),
               constant_matches_formula=bool(np.allclose(const[0], PB.tsp_constant(W, A))),
               field_ok=bool(np.allclose(pr.field(S[:64]), S[:64] @ J + b[None, :])),
               optimal_tour=int(Lopt), feasible_min_energy_is_4Lopt=bool(np.isclose(emin_feas, 4 * Lopt + const[0])),
               no_infeasible_below_feasible_min=bool(np.all(e[~feas] >= emin_feas)),
               ground_states_feasible=bool(np.all(feas[gs])),
               decode_ok=bool(np.all(pr.evaluate(S[feas][:200])['feasible'])))
    if A_rule > 1:
        res['ok'] = bool(all(v for k, v in res.items() if isinstance(v, bool)))
    else:   # A = max W (MQC): degenerate infeasible ground states allowed, none strictly below the feasible minimum
        res['ok'] = bool(all(v for k, v in res.items() if isinstance(v, bool) and k != 'ground_states_feasible'))
    return res


def part_A():
    out = []
    for n, p, P, seed in ((10, 0.4, 4, 1), (12, 0.35, 4, 2), (12, 0.35, 2, 3), (14, 0.3, 8, 4), (14, 0.3, 4, 5)):
        r = brute_gpp(n, p, P, seed); out.append(r)
        log('A gpp', json.dumps(r, default=str))
    for n, seed, ar in ((3, 1, 1.0), (4, 2, 1.0), (4, 3, 1.5), (4, 4, 1.0)):
        r = brute_tsp(n, seed, ar); out.append(r)
        log('A tsp', json.dumps(r, default=str))
    OUT['A_formulations'] = out
    return all(r['ok'] for r in out)


# ============================================================ B. bit-identity with the G-set kernels (b = 0, gamma = 0)
def part_B():
    import engine as E
    import others as O
    import synth
    import run_gset as RG
    bkv = json.loads((RES / 'gset_20261007' / 'bkv.json').read_text())['instances']
    cases = []; ok = True
    gs = [(E.Graph('G1'), 1), (E.Graph('G14'), 14), (synth.random_graph(300, 2700, True, 1), None),
          (synth.planar_like(300, False, 5), None)]
    for g, gi in gs:
        P = PB.mcp(gi, bkv[f'G{gi}']) if gi else None
        if P is None:
            rows = np.concatenate([g.ei, g.ej]); cols = np.concatenate([g.ej, g.ei]); vals = -np.concatenate([g.w, g.w])
            P = PB.Problem(g.name, 'mcp', g.N, rows, cols, vals, 0.0, np.zeros(g.N), g.sigma,
                           dict(ei=g.ei, ej=g.ej, w=g.w, ref=1, target=1))
        assert np.array_equal(P.indptr, g.indptr) and np.array_equal(P.indices, g.indices) and np.array_equal(P.data, g.data)
        a = g.sigma; B = 8
        cfgs = [dict(family='plain', q=1.44 * a, T0=1.35 * a, tfin=0.672 * a, S=120),
                dict(family='tec', q=1.44 * a, jv=-0.18 * a, T0=1.35 * a, tfin=0.672 * a, S=120),
                dict(family='tecT', q=1.44 * a, kappa=0.5, ramp=True, T0=1.35 * a, tfin=0.672 * a, S=120),
                dict(family='onsager', q=1.44 * a, lam=0.5 * g.lam_factor, ramp=True, T0=1.35 * a, tfin=0.672 * a, S=120),
                dict(family='apc', q_reset=2.88 * a, r_q=0.9, q_lim=1.44 * a, T0=1.35 * a, tfin=0.672 * a, S=120),
                dict(family='SA', T0=0.671 * a, T1=0.0447 * a, S=60),
                dict(family='ReAIM', kset=(26, 51, 102, 205), T1=0.05, S=100),
                dict(family='dSB', dt=1.0, xi=1.0, S=120), dict(family='bSB', dt=0.75, xi=1.0, S=120),
                dict(family='aSB', dt=0.5, xi=1.0, S=120)]
        for cfg in cfgs:
            seed = np.random.SeedSequence([1, 2, 3])
            s1, f1, fin1 = SV.run_method(P, cfg, B, seed)
            seed = np.random.SeedSequence([1, 2, 3])
            rng = np.random.default_rng(seed); fam = cfg['family']
            if fam in SV.ENGINE_FAMILIES:
                s0 = rng.choice(np.array([-1.0, 1.0], dtype=np.float32), size=(B, g.N)); s2, f2 = E.run(g, s0, cfg, rng, cfg['tfin'])
            elif fam == 'SA':
                s2 = O.sa(g, B, cfg, rng); f2 = None
            elif fam == 'ReAIM':
                s2 = O.reaim(g, B, cfg, rng); f2 = None
            elif fam in ('bSB', 'dSB'):
                s2 = O.sb(g, B, cfg, rng); f2 = None
            else:
                s2, fin2 = O.asb(g, B, cfg, rng); f2 = None
            same = bool(np.array_equal(s1, s2) and (f2 is None or np.array_equal(f1, f2)))
            cut_same = bool(np.array_equal(P.evaluate(s1)['value'], g.cut(s2)))
            ok &= same and cut_same
            cases.append(dict(graph=g.name, family=fam, identical=same, cut_identical=cut_same))
            log(f'B {g.name:24s} {fam:8s} identical={same} cut={cut_same}')
    # run_gset.run_method itself on G2 (an MCP-extension instance) for one grid point of each other method
    G2 = RG.graph(2)
    for meth in RG.OTHERS:
        cfg = RG.grid(meth, G2, 2, 250)[5]
        c1, _, _ = RG.run_method(G2, cfg, 4, np.random.SeedSequence([9, 9]))
        P2 = PB.mcp(2, bkv['G2'])
        s, _, _ = SV.run_method(P2, cfg, 4, np.random.SeedSequence([9, 9]))
        same = bool(np.array_equal(P2.evaluate(s)['value'], c1))
        ok &= same
        cases.append(dict(graph='G2', family=meth, via='run_gset.run_method', identical=same))
        log(f'B G2 via run_gset {meth:10s} identical={same}')
    OUT['B_gset_identity'] = cases
    return ok


# ============================================================ C. bias and uniform coupling against dense references
def abl_run_bias(J, b, s0, cfg, rng, tfin):
    """research/ablation_20261005/abl.py::run, copied, with the schedule's final temperature as a parameter and the single
    change h = s @ J + b at initialization (marked)."""
    fam = cfg['family']; S = cfg['S']; T = SV.sched(cfg['T0'], S, tfin)
    s = s0.astype(np.float32).copy(); h = s @ J + b[None, :].astype(np.float32)          # <- bias
    B, N = s.shape
    flips = np.zeros(B); s_prev = None; c_prev = None; T_prev = None
    q = np.full((B, N), float(cfg.get('q', 0.0)), np.float32)
    if fam == 'apc':
        q[:] = cfg['q_reset']
    for t, Tt in enumerate(T):
        field = h
        if s_prev is not None:
            if fam == 'onsager':
                field = h - cfg['lam'] * SV.ramp_factor(t, S, cfg['ramp']) * c_prev * s_prev
            elif fam == 'tecT':
                field = h - cfg['kappa'] * T_prev * SV.ramp_factor(t, S, cfg['ramp']) * s_prev
            elif fam == 'tec':
                field = h + cfg['jv'] * s_prev
        z = s * field + q
        flip = np.clip(z / (4 * Tt) + 0.5, 0, 1) < rng.random(s.shape, dtype=np.float32)
        if fam in ('onsager',):
            c_prev = ((z > -2 * Tt) & (z < 2 * Tt)).sum(1, keepdims=True) / (2 * Tt)
        if fam == 'apc':
            q = np.where(flip, np.float32(cfg['q_reset']), np.maximum(q * np.float32(cfg['r_q']), np.float32(cfg['q_lim'])))
        d = np.where(flip, -2 * s, 0).astype(np.float32)
        s_prev = s.copy(); T_prev = Tt
        s += d; h += d @ J; flips += flip.sum(1)
    return s, flips


def sa_dense_ref(J, b, s0, Tsched, seeds):
    """Pure numpy SA with the legacy MT19937 stream (np.random.seed / np.random.random), dense J, bias."""
    B, N = s0.shape; s = s0.astype(np.float64).copy()
    for r in range(B):
        np.random.seed(int(seeds[r])); sb = s[r].copy(); h = J @ sb + b
        for Tt in Tsched:
            for i in range(N):
                dE = 2.0 * sb[i] * h[i]
                if dE <= 0.0 or np.random.random() < math.exp(-dE / Tt):
                    sb[i] = -sb[i]; h += 2.0 * sb[i] * J[i]
        s[r] = sb
    return s


def dense_problem(P):
    """Same couplings as P with the uniform part folded into explicit CSR entries (gamma = 0)."""
    J = P.dense(); r, c = np.nonzero(J)
    return PB.Problem(P.name + '_dense', P.kind, P.N, r, c, J[r, c], 0.0, P.b.astype(np.float64), P.sigma_T, P.meta)


def part_C():
    cases = []; ok = True
    rng0 = np.random.default_rng(11)
    n = 60; E = [(i, j) for i in range(n) for j in range(i + 1, n) if rng0.random() < 0.15]
    ei = np.array([e[0] for e in E]); ej = np.array([e[1] for e in E])
    gp = PB.gpp_from_edges('c_gpp60', n, ei, ej, np.ones(len(E), np.int64), 4, ref=1, target=1)
    pts = rng0.random((6, 2)) * 100
    W = np.rint(np.hypot(*(pts[:, None, :] - pts[None, :, :]).transpose(2, 0, 1))).astype(np.int64); np.fill_diagonal(W, 0)
    W[(W == 0) & ~np.eye(6, dtype=bool)] = 1
    tp = PB.tsp_from_matrix('c_tsp6', W, float(W.max()), ref=1, target=1)
    for P in (gp, tp):
        J = P.dense().astype(np.float32); b = P.b
        a = P.sigma_T; B = 6
        Pd = dense_problem(P)
        cfgs = [dict(family='plain', q=1.0 * a, T0=1.5 * a, tfin=0.3 * a, S=80),
                dict(family='tec', q=1.0 * a, jv=-0.2 * a, T0=1.5 * a, tfin=0.3 * a, S=80),
                dict(family='tecT', q=1.0 * a, kappa=0.01, ramp=True, T0=1.5 * a, tfin=0.3 * a, S=80),
                dict(family='onsager', q=1.0 * a, lam=0.05 * P.stats()['lam_factor'], ramp=True, T0=1.5 * a, tfin=0.3 * a, S=80),
                dict(family='apc', q_reset=2.0 * a, r_q=0.9, q_lim=0.5 * a, T0=1.5 * a, tfin=0.3 * a, S=80)]
        for cfg in cfgs:
            r1 = np.random.default_rng(5); s0 = r1.choice(np.array([-1.0, 1.0], dtype=np.float32), size=(B, P.N))
            s1, f1 = SV.engine_run(P, s0, cfg, r1, cfg['tfin'])
            r2 = np.random.default_rng(5); s0b = r2.choice(np.array([-1.0, 1.0], dtype=np.float32), size=(B, P.N))
            s2, f2 = abl_run_bias(J, b, s0b, cfg, r2, cfg['tfin'])
            same = bool(np.array_equal(s1, s2) and np.array_equal(f1, f2))
            ok &= same; cases.append(dict(problem=P.name, family=cfg['family'], ref='abl.run+bias (dense)', identical=same,
                                          mean_flips=float(f1.mean())))
            log(f'C {P.name:8s} {cfg["family"]:8s} vs dense abl.run+bias identical={same} flips={f1.mean():.0f}')
        # SA vs pure numpy legacy-MT reference
        cfg = dict(family='SA', T0=1.0 * a, T1=0.05 * a, S=30)
        r1 = np.random.default_rng(6); s_a = SV.sa(P, 3, cfg, r1)
        r2 = np.random.default_rng(6); s0 = r2.choice(np.array([-1.0, 1.0], dtype=np.float32), size=(3, P.N))
        seeds = r2.integers(0, 2 ** 31 - 1, size=3)
        Ts = cfg['T0'] * (cfg['T1'] / cfg['T0']) ** (np.arange(cfg['S']) / (cfg['S'] - 1))
        s_b = sa_dense_ref(P.dense(), P.b.astype(np.float64), s0, Ts, seeds)
        same = bool(np.array_equal(s_a, s_b.astype(np.float32))); ok &= same
        cases.append(dict(problem=P.name, family='SA', ref='numpy dense legacy-MT', identical=same))
        log(f'C {P.name:8s} SA vs numpy dense reference identical={same}')
        # ReAIM, dSB, bSB, aSB, engine: uniform part vs folded explicit entries (and bias kept)
        for cfg in (dict(family='ReAIM', kset=(2, 4, 8, 16), T1=0.05, S=60), dict(family='dSB', dt=1.0, xi=1.0, S=80),
                    dict(family='SA', T0=1.0 * a, T1=0.05 * a, S=30), cfgs[3]):
            sA, _, _ = SV.run_method(P, cfg, 4, np.random.SeedSequence([4, 4]))
            sB, _, _ = SV.run_method(Pd, cfg, 4, np.random.SeedSequence([4, 4]))
            same = bool(np.array_equal(sA, sB)); ok &= same
            cases.append(dict(problem=P.name, family=cfg['family'], ref='gamma folded into CSR', identical=same))
            log(f'C {P.name:8s} {cfg["family"]:8s} gamma vs folded identical={same}')
        # ReAIM with bias vs dense numpy Batch (same code path with dense products)
        cfg = dict(family='ReAIM', kset=(2, 4, 8, 16), T1=0.05, S=60)
        sA, _, _ = SV.run_method(P, cfg, 3, np.random.SeedSequence([8]))
        sR = reaim_dense_ref(P, cfg, 3, np.random.SeedSequence([8]))
        same = bool(np.array_equal(sA, sR)); ok &= same
        cases.append(dict(problem=P.name, family='ReAIM', ref='dense numpy Batch with bias', identical=same))
        log(f'C {P.name:8s} ReAIM vs dense reference identical={same}')
        # dSB with bias vs dense float64 reference (integer sums: exact)
        cfg = dict(family='dSB', dt=1.0, xi=1.0, S=60)
        sA, _, _ = SV.run_method(P, cfg, 3, np.random.SeedSequence([12]))
        sR = dsb_dense_ref(P, cfg, 3, np.random.SeedSequence([12]))
        same = bool(np.array_equal(sA, sR)); ok &= same
        cases.append(dict(problem=P.name, family='dSB', ref='dense reference with bias', identical=same))
        log(f'C {P.name:8s} dSB vs dense reference identical={same}')
    OUT['C_bias_gamma'] = cases
    return ok


def reaim_dense_ref(P, cfg, B, seed):
    rng = np.random.default_rng(seed)
    J = P.dense().astype(np.float32); b = P.b.astype(np.float32)
    S = cfg['S']; kset = cfg['kset']; t0, t1 = 1.0, cfg['T1']; it_trial, it_run = 32, 96; N = P.N; R = len(kset)
    alpha = (t1 / t0) ** (1.0 / (S - 1))

    class Bt:
        def __init__(s_, x, fifo, h=None):
            s_.x, s_.fifo = x, fifo; s_.h = (x @ J + b[None, :]) if h is None else h

        def energy(s_):
            return -0.5 * np.einsum('bi,bi->b', s_.x, s_.h + b[None, :]) if np.any(b) else -0.5 * np.einsum('bi,bi->b', s_.x, s_.h)

        def step(s_, k, T):
            x, h = s_.x, s_.h
            d = 2.0 * x * h; nN = (d < 0).sum(1); q = np.max(s_.fifo, axis=1)
            pos = np.where(d > 0, d, np.inf).min(1); dmin, dmax = d.min(1), d.max(1)
            th_pm = np.where(np.isfinite(pos), pos, dmax); th_mn = dmin + rng.random(len(x)) * T * (dmax - dmin)
            th = np.where(nN > q, th_pm, th_mn); C = d <= th[:, None]
            keys = np.where(C, rng.random(C.shape), -1.0)
            kk = np.broadcast_to(np.asarray(k), (len(x),)); kmax = int(kk.max())
            idx = np.argpartition(-keys, kmax - 1, axis=1)[:, :kmax]
            sel = np.take_along_axis(keys, idx, 1) >= 0; sel &= np.arange(kmax)[None, :] < kk[:, None]
            rows = np.nonzero(sel); flat = idx[rows]
            delta = np.zeros_like(x); delta[rows[0], flat] = -2.0 * x[rows[0], flat]
            x += delta; s_.h = h + delta @ J
            s_.fifo = np.concatenate([s_.fifo[:, 1:], nN[:, None]], axis=1)
    x = rng.choice(np.array([-1.0, 1.0], dtype=np.float32), size=(B, N))
    cur = Bt(x, np.full((B, 20), N, dtype=np.int64)); t, T = 0, t0
    while t < S:
        n = min(it_trial, S - t)
        rep = Bt(np.repeat(cur.x, R, 0), np.repeat(cur.fifo, R, 0), h=np.repeat(cur.h, R, 0))
        kvec = np.tile(np.array(kset), B)
        for _ in range(n):
            rep.step(kvec, T); T *= alpha
        pick = rep.energy().reshape(B, R).argmin(1); t += n
        sel = np.arange(B) * R + pick
        cur = Bt(rep.x[sel].copy(), rep.fifo[sel].copy(), h=rep.h[sel].copy()); kbest = np.array(kset)[pick]
        n = min(it_run, S - t)
        for _ in range(n):
            cur.step(kbest, T); T *= alpha
        t += n
    return cur.x


def dsb_dense_ref(P, cfg, B, seed):
    rng = np.random.default_rng(seed)
    J = P.dense(); b = P.b.astype(np.float64)
    S = cfg['S']; dt = np.float32(cfg['dt'])
    c0 = np.float32(cfg['xi'] * 0.5 / (P.stats()['sigmaJ_goto'] * math.sqrt(P.N)))
    x = (rng.random((B, P.N), dtype=np.float32) - 0.5) * 0.2
    y = (rng.random((B, P.N), dtype=np.float32) - 0.5) * 0.2
    for t in range(S):
        ka = np.float32(-(1.0 - 1.0 * t / S))
        xs = np.where(x >= 0, 1.0, -1.0)
        jx = (xs @ J + b[None, :]).astype(np.float32)        # exact integers
        y = (y + dt * (ka * x + c0 * jx)).astype(np.float32)
        xv = (x + np.float32(dt * 1.0) * y).astype(np.float32)
        out = np.abs(xv) > 1
        x = np.where(out, np.where(xv > 0, np.float32(1), np.float32(-1)), xv).astype(np.float32)
        y = np.where(out, np.float32(0), y).astype(np.float32)
    return np.where(x >= 0, 1.0, -1.0).astype(np.float32)


def main():
    t0 = time.time()
    okA = part_A(); okC = part_C(); okB = part_B()
    OUT.update(ok_A=okA, ok_B=okB, ok_C=okC, all_ok=bool(okA and okB and okC), sec=time.time() - t0)
    import platform, numba
    OUT['platform'] = dict(machine=platform.machine(), python=platform.python_version(), numpy=np.__version__,
                           numba=numba.__version__, node=platform.node())
    tag = os.environ.get('VERIFY_TAG', platform.node().split('.')[0])
    (HERE / f'verify_{tag}.json').write_text(json.dumps(OUT, indent=1, default=str) + '\n')
    log('A', okA, 'B', okB, 'C', okC, '->', 'ALL OK' if OUT['all_ok'] else 'FAIL')
    return 0 if OUT['all_ok'] else 1


if __name__ == '__main__':
    sys.exit(main())
