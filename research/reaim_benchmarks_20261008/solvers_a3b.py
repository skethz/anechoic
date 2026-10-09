"""Addendum 3.1: the output rules the papers describe.
  * APC-SCA (Okonogi et al. 2023, Algorithms 1-2): "Output: optimized spin states: argmin_{1<=s<=S+1} H(sigma(s))",
    the lowest-energy state visited, including the initial state.            -> engine_sig_best
  * ReAIM ASA (Chiang et al. 2024, Algorithm 3): "Output: best spin state x_best, best Hamiltonian H_best", updated after
    every run phase (lines 6-8).                                               -> reaim_best
Both perform exactly the dynamics (and RNG use) of solvers_a3.engine_sig / solvers_a2.reaim_a2; only the returned
state differs (verify_a3b.py). Every other method returns its final state, as its paper describes."""
import math

import numpy as np
import numba as nb

import solvers as SV
import solvers_a3 as SV3


@nb.njit(cache=True)
def _energy(s, h, Mv, gamma, b):
    """E = -1/2 sum_i s_i F_i - 1/2 sum_i b_i s_i, F_i = h_i + gamma (M - s_i) (h includes b): -sum_{i<j} J s s - sum b s."""
    B, N = s.shape
    E = np.zeros(B)
    for r in range(B):
        acc = 0.0
        for i in range(N):
            F = np.float64(h[r, i]) + gamma * (Mv[r] - s[r, i])
            acc += s[r, i] * F + b[i] * s[r, i]
        E[r] = -0.5 * acc
    return E


def engine_sig_best(P, s0, cfg, rng):
    """solvers_a3.engine_sig with the APC output rule; returns (best spins, flips, final spins)."""
    S = cfg['S']; T = SV.sched(cfg['T0'], S, cfg['tfin'])
    s = s0.astype(np.float64).copy()
    h = SV.field_csr_b(s.astype(np.float32), P.indptr, P.indices, P.data, P.b)
    B, N = s.shape
    Mv = s.sum(1)
    flips = np.zeros(B); sp = s.copy()
    apc = cfg['family'] == 'apc_sig'; tec = cfg['family'] == 'tec_sig'
    q = np.full((B, N), float(cfg['q_reset']) if apc else float(cfg.get('q', 0.0)))
    flipbuf = np.zeros(N, np.int64)
    b64 = P.b.astype(np.float64)
    bestE = _energy(s, h, Mv, float(P.gamma), b64); best = s.copy()
    for t, Tt in enumerate(T):
        u = rng.random(s.shape)
        jv = float(cfg['jv']) if (tec and t > 0) else 0.0
        SV3._step_sig(apc, tec and t > 0, s, h, Mv, float(P.gamma), sp, q, u, float(Tt), float(cfg['cb']), jv, flips,
                      P.indptr, P.indices, P.data, float(cfg.get('q_reset', 0.0)), float(cfg.get('r_q', 0.0)),
                      float(cfg.get('q_lim', 0.0)), flipbuf)
        E = _energy(s, h, Mv, float(P.gamma), b64)
        imp = E < bestE
        if imp.any():
            best[imp] = s[imp]; bestE[imp] = E[imp]
    return best.astype(np.float32), flips, s.astype(np.float32)


def reaim_best(P, B, cfg, rng):
    """solvers_a2.reaim_a2 with Algorithm 3's output: the lowest-energy state among the run-phase end states."""
    F = {'max': np.max, 'min': np.min}[cfg.get('F', 'max')]
    S = cfg['S']; kset = cfg['kset']; t0, t1 = float(cfg.get('T0', 1.0)), cfg['T1']; it_trial, it_run = 32, 96
    N = P.N; R = len(kset)
    alpha = (t1 / t0) ** (1.0 / (S - 1))
    x = rng.choice(np.array([-1.0, 1.0], dtype=np.float32), size=(B, N))
    cur = SV.Batch(P, x, np.full((B, 20), N, dtype=np.int64), F)
    t, T = 0, t0
    Hbest = np.full(B, np.inf); xbest = cur.x.copy()
    while t < S:
        n = min(it_trial, S - t)
        rep = SV.Batch(P, np.repeat(cur.x, R, 0), np.repeat(cur.fifo, R, 0), F, h=np.repeat(cur.h, R, 0),
                       M=np.repeat(cur.M, R, 0))
        kvec = np.tile(np.array(kset), B)
        for u in range(n):
            rep.step(kvec, T, rng); T *= alpha
        Elast = rep.energy().reshape(B, R)
        pick = Elast.argmin(1)
        t += n
        sel = np.arange(B) * R + pick
        cur = SV.Batch(P, rep.x[sel].copy(), rep.fifo[sel].copy(), F, h=rep.h[sel].copy(), M=rep.M[sel].copy())
        kbest = np.array(kset)[pick]
        n = min(it_run, S - t)
        for u in range(n):
            cur.step(kbest, T, rng); T *= alpha
        t += n
        H = cur.energy()
        imp = H < Hbest
        if imp.any():
            xbest[imp] = cur.x[imp]; Hbest[imp] = H[imp]
    return xbest, cur.x


def run_method(P, cfg, B, seed):
    """APC-SCA and ReAIM with their papers' output rules; everything else as solvers_a3.run_method (final state)."""
    if cfg['family'] == 'apc_sig':
        rng = np.random.default_rng(seed)
        s0 = rng.choice(np.array([-1.0, 1.0], dtype=np.float32), size=(B, P.N))
        best, flips, final = engine_sig_best(P, s0, cfg, rng)
        return best, flips, None
    if cfg['family'] == 'ReAIM':
        rng = np.random.default_rng(seed)
        best, final = reaim_best(P, B, cfg, rng)
        return best, None, None
    return SV3.run_method(P, cfg, B, seed)
