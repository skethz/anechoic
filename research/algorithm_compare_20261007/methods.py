"""Algorithm-level comparison on K2000: Ising energy H(s) = -1/2 s^T J s versus Monte Carlo step, for every method.

Step convention (each step is one full field update, O(N^2) work):
  SA       one sweep of N sequential single-spin Metropolis updates (spins 0..N-1 in order)
  SCA      one synchronous update of all N spins (STATICA's clipped probability, self-coupling penalty q)
           families: plain, tec (+J_v s(t-1)), apc (adaptive penalty), onsager (ours, -lam c(t-1) s(t-1)), tecT (ours, -kappa T s(t-1))
  bSB/dSB  one time step of ballistic / discrete simulated bifurcation (Goto et al., Sci. Adv. 2021): one J x product
Energy: H = sumw - 2 cut, so cut >= 33,000 <=> H <= -67,040 and the best-known cut 33,337 <=> H = -67,714.
All methods start from independent uniformly random spins (SB from small random amplitudes) with per-run seeds.
"""
import math
import sys
from pathlib import Path

import numpy as np
import numba as nb

sys.path.insert(0, str(__import__('pathlib').Path(__file__).resolve().parent / '../ablation_20261005'))
sys.path.insert(0, str(__import__('pathlib').Path(__file__).resolve().parent / '../reaim_reproduction_20261003'))
import abl  # noqa: E402  (clipped parallel SCA families, the same model used for every ablation and hardware schedule)
import pilot as reaim  # noqa: E402  (noise-free ReAIM ASA, Algorithms 2-3, as reproduced on 3 October 2026)

m = abl.m


def energy_from_field(s, h):
    return -0.5 * np.einsum('bi,bi->b', s, h)


# ---------------------------------------------------------------- SCA families with per-step energy trajectory
def sca_traj(J, s0, cfg, rng):
    """abl.run with the per-step energy recorded (same arithmetic as abl.run)."""
    fam = cfg['family']; S = cfg['S']; T = m.sched(cfg['T0'], S)
    s = s0.astype(np.float32).copy(); h = s @ J
    B, N = s.shape
    E = np.empty((S + 1, B), np.float64); E[0] = energy_from_field(s, h)
    s_prev = None; c_prev = None; T_prev = None
    q = np.full((B, N), float(cfg.get('q', 0.0)), np.float32)
    if fam == 'apc':
        q[:] = cfg['q_reset']
    for t, Tt in enumerate(T):
        field = h
        if s_prev is not None:
            if fam == 'onsager':
                field = h - cfg['lam'] * abl.ramp_factor(t, S, cfg['ramp']) * c_prev * s_prev
            elif fam == 'tecT':
                field = h - cfg['kappa'] * T_prev * abl.ramp_factor(t, S, cfg['ramp']) * s_prev
            elif fam == 'tec':
                field = h + cfg['jv'] * s_prev
        z = s * field + q
        flip = np.clip(z / (4 * Tt) + 0.5, 0, 1) < rng.random(s.shape, dtype=np.float32)
        if fam == 'onsager':
            c_prev = ((z > -2 * Tt) & (z < 2 * Tt)).sum(1, keepdims=True) / (2 * Tt)
        if fam == 'apc':
            q = np.where(flip, np.float32(cfg['q_reset']), np.maximum(q * np.float32(cfg['r_q']), np.float32(cfg['q_lim'])))
        d = np.where(flip, -2 * s, 0).astype(np.float32)
        s_prev = s.copy(); T_prev = Tt
        s += d; h += d @ J
        E[t + 1] = energy_from_field(s, h)
    return s, E


# ---------------------------------------------------------------- SA (sequential Metropolis sweeps), numba, parallel over runs
@nb.njit(parallel=True, fastmath=False, cache=True)
def _sa_kernel(J, s, Tsched, seeds, E):
    B, N = s.shape
    S = Tsched.shape[0]
    for b in nb.prange(B):
        np.random.seed(seeds[b])
        sb = s[b].copy()
        h = np.zeros(N)
        for i in range(N):
            acc = 0.0
            for j in range(N):
                acc += J[i, j] * sb[j]
            h[i] = acc
        e = 0.0
        for i in range(N):
            e += sb[i] * h[i]
        E[0, b] = -0.5 * e
        for t in range(S):
            Tt = Tsched[t]
            for i in range(N):
                dE = 2.0 * sb[i] * h[i]
                if dE <= 0.0 or np.random.random() < math.exp(-dE / Tt):
                    sb[i] = -sb[i]
                    d = 2.0 * sb[i]
                    for j in range(N):
                        h[j] += d * J[i, j]   # J symmetric: row access
            e = 0.0
            for i in range(N):
                e += sb[i] * h[i]
            E[t + 1, b] = -0.5 * e
        s[b] = sb


def sa_traj(J, s0, cfg, seeds):
    S = cfg['S']
    Tsched = cfg['T0'] * (cfg['T1'] / cfg['T0']) ** (np.arange(S) / max(1, S - 1))
    s = s0.astype(np.float64).copy(); E = np.zeros((S + 1, s.shape[0]))
    _sa_kernel(J.astype(np.float64), s, Tsched.astype(np.float64), np.asarray(seeds, np.int64), E)
    return s.astype(np.float32), E


# ---------------------------------------------------------------- simulated bifurcation (Goto et al. 2021)
def sb_traj(J, B, cfg, rng):
    """Ballistic (bSB) or discrete (dSB) SB: y += dt[-(a0 - a(t)) x + c0 J x'], x += dt a0 y, inelastic walls |x| <= 1.
    x' = x (bSB) or sign(x) (dSB); a(t) rises linearly from 0 to a0 over S steps; c0 = xi0 * 0.5 / (sigma sqrt(N))."""
    S = cfg['S']; dt = cfg['dt']; a0 = 1.0
    N = J.shape[0]
    sigma = float(np.sqrt((J ** 2).sum() / (N * (N - 1))))
    c0 = cfg.get('xi', 1.0) * 0.5 / (sigma * math.sqrt(N))
    x = (rng.random((B, N), dtype=np.float32) - 0.5) * 0.2
    y = (rng.random((B, N), dtype=np.float32) - 0.5) * 0.2
    E = np.empty((S + 1, B))
    s = np.where(x >= 0, 1.0, -1.0).astype(np.float32); E[0] = energy_from_field(s, s @ J)
    disc = cfg['family'] == 'dSB'
    for t in range(S):
        at = a0 * t / S
        xs = np.where(x >= 0, 1.0, -1.0).astype(np.float32) if disc else x
        y += dt * (-(a0 - at) * x + c0 * (xs @ J))
        x += dt * a0 * y
        wall = np.abs(x) > 1
        x = np.where(wall, np.sign(x), x).astype(np.float32); y = np.where(wall, 0.0, y).astype(np.float32)
        s = np.where(x >= 0, 1.0, -1.0).astype(np.float32)
        E[t + 1] = energy_from_field(s, s @ J)
    return s, E


def asb_traj(J, B, cfg, rng):
    """Adiabatic SB (Goto et al. 2019), same integrator as gpu/prior_compare_20260929/sb_gpu.py: M Kerr sub-steps, then
    y += dt*xi*J x; pump p = t/S; x(0) = 0, y(0) = 0.1*u*s0; xi = xi_mult*0.7/sqrt(N)."""
    S = cfg['S']; dt = cfg['dt']; M = cfg.get('M', 5); N = J.shape[0]
    xi = np.float32(cfg['xi'] * 0.7 / math.sqrt(N))
    s0 = rng.choice(np.array([-1.0, 1.0], dtype=np.float32), size=(B, N))
    x = np.zeros((B, N), np.float32); y = (s0 * 0.1 * rng.random((B, N), dtype=np.float32)).astype(np.float32)
    E = np.empty((S + 1, B)); E[0] = energy_from_field(s0, s0 @ J)
    delta = np.float32(dt / M); finite = True; s = s0
    for t in range(1, S + 1):
        p = np.float32(t / S)
        for _ in range(M):
            x = x + y * delta
            y = y - ((x * x + 1 - p) * x) * delta
        y = y + (x @ J) * (xi * np.float32(dt))
        if not (np.isfinite(x).all() and np.isfinite(y).all()):
            finite = False; x = np.nan_to_num(x); y = np.nan_to_num(y)
        s = np.where(x >= 0, 1.0, -1.0).astype(np.float32)
        E[t] = energy_from_field(s, s @ J)
    return s, E, finite


def reaim_traj(J, B, cfg, rng):
    """pilot.asa() with the energy of the search's own trajectory recorded at every iteration. During a trial phase the
    trajectory is the replica that the phase then selects (lowest energy at the end of the phase), as in Algorithm 3."""
    S = cfg['S']; kset = cfg['kset']; t0, t1 = 1.0, cfg['T1']; it_trial, it_run = 32, 96
    N = J.shape[0]; R = len(kset)
    alpha = (t1 / t0) ** (1.0 / (S - 1))
    x = rng.choice(np.array([-1.0, 1.0], dtype=np.float32), size=(B, N))
    cur = reaim.Batch(J, x, np.full((B, 20), N, dtype=np.int64), np.max)
    E = np.empty((S + 1, B)); E[0] = cur.energy()
    t, T = 0, t0
    while t < S:
        n = min(it_trial, S - t)
        rep = reaim.Batch(J, np.repeat(cur.x, R, 0), np.repeat(cur.fifo, R, 0), np.max)
        kvec = np.tile(np.array(kset), B)
        Etr = np.empty((n, B, R))
        for u in range(n):
            rep.step(kvec, T, rng); T *= alpha
            Etr[u] = rep.energy().reshape(B, R)
        pick = Etr[-1].argmin(1)
        E[t + 1:t + 1 + n] = Etr[:, np.arange(B), pick]
        t += n
        sel = np.arange(B) * R + pick
        cur = reaim.Batch(J, rep.x[sel].copy(), rep.fifo[sel].copy(), np.max)
        kbest = np.array(kset)[pick]
        n = min(it_run, S - t)
        for u in range(n):
            cur.step(kbest, T, rng); T *= alpha
            E[t + 1 + u] = cur.energy()
        t += n
    return cur.x, E


def run_method(J, sumw, cfg, B, seed):
    rng = np.random.default_rng(seed)
    fam = cfg['family']
    finite = True
    if fam in ('bSB', 'dSB'):
        s, E = sb_traj(J, B, cfg, rng)
    elif fam == 'aSB':
        s, E, finite = asb_traj(J, B, cfg, rng)
    elif fam == 'ReAIM':
        s, E = reaim_traj(J, B, cfg, rng)
    else:
        s0 = rng.choice(np.array([-1.0, 1.0], dtype=np.float32), size=(B, len(J)))
        if fam == 'SA':
            s, E = sa_traj(J, s0, cfg, rng.integers(0, 2 ** 31 - 1, size=B))
        else:
            s, E = sca_traj(J, s0, cfg, rng)
    cuts = (sumw - E[-1]) / 2.0
    if not finite:
        cuts = np.full_like(cuts, -np.inf)   # non-finite dynamics: configuration invalid
    return E, cuts
