"""The non-engine methods of research/algorithm_compare_20261007/methods.py (not edited), adapted to sparse couplings.

  SA      methods._sa_kernel with CSR rows (adding J_ij = 0 terms is a no-op), same per-run np.random seeding and draws:
          decisions are bit-identical (verify_others.py).
  ReAIM   pilot.Batch / methods.reaim_traj with the field update as an exact integer scatter: bit-identical.
  dSB     methods.sb_traj with sign(x) couplings: the J.sign(x) products are exact integers, float32 elsewhere in the
          same operation order: bit-identical up to float32 rounding of identical expressions (verified).
  bSB     same integrator; J.x for real x is summed over CSR rows in float32, so the summation order (and hence the last
          bits) differs from dense BLAS: algorithmically identical, not bit-identical.
  aSB     methods.asb_traj (Goto et al. 2019 integrator, M = 5 Kerr sub-steps) with xi0 = xi * 0.7 / (sigma_J sqrt(N)),
          sigma_J = sqrt(sum J^2 / (N(N-1))): the published normalisation. methods.py hard-wires sigma_J = 1 (exact
          for K2000), which would make the coupling ~sqrt(N/d) times too weak on sparse graphs.
All return final spins; cut is computed exactly by Graph.cut.
"""
import math

import numpy as np
import numba as nb

import engine as E


# ------------------------------------------------------------------ SA (sequential Metropolis sweeps)
@nb.njit(cache=True)
def _sa_kernel_csr(indptr, indices, data, s, Tsched, seeds):
    B, N = s.shape
    S = Tsched.shape[0]
    for b in range(B):
        np.random.seed(seeds[b])
        sb = s[b].copy()
        h = np.zeros(N)
        for i in range(N):
            acc = 0.0
            for k in range(indptr[i], indptr[i + 1]):
                acc += data[k] * sb[indices[k]]
            h[i] = acc
        for t in range(S):
            Tt = Tsched[t]
            for i in range(N):
                dE = 2.0 * sb[i] * h[i]
                if dE <= 0.0 or np.random.random() < math.exp(-dE / Tt):
                    sb[i] = -sb[i]
                    d = 2.0 * sb[i]
                    for k in range(indptr[i], indptr[i + 1]):
                        h[indices[k]] += d * data[k]
        s[b] = sb


def sa(g, B, cfg, rng):
    """Same RNG use as methods.run_method for SA: s0 from rng.choice, then per-run seeds from rng.integers."""
    s0 = rng.choice(np.array([-1.0, 1.0], dtype=np.float32), size=(B, g.N))
    seeds = rng.integers(0, 2 ** 31 - 1, size=B)
    S = cfg['S']
    Tsched = cfg['T0'] * (cfg['T1'] / cfg['T0']) ** (np.arange(S) / max(1, S - 1))
    s = s0.astype(np.float64).copy()
    _sa_kernel_csr(g.indptr, g.indices, g.data.astype(np.float64), s, Tsched.astype(np.float64), np.asarray(seeds, np.int64))
    return s.astype(np.float32)


# ------------------------------------------------------------------ simulated bifurcation (bSB / dSB), Goto et al. 2021
def sigma_J(g):
    return math.sqrt(2.0 * g.m / (g.N * (g.N - 1)))     # sqrt(sum_ij J_ij^2 / (N(N-1))), unit weights


@nb.njit(cache=True)
def _sb_kernel(indptr, indices, data, x, y, S, dt, c0, disc):
    B, N = x.shape
    a0 = 1.0
    jx = np.zeros(N, np.float32); xs = np.zeros(N, np.float32)
    dtf = np.float32(dt); c0f = np.float32(c0); dta = np.float32(dt * a0)
    one = np.float32(1.0); zero = np.float32(0.0)
    for b in range(B):
        for t in range(S):
            at = a0 * t / S
            ka = np.float32(-(a0 - at))
            for i in range(N):
                if disc:
                    xs[i] = one if x[b, i] >= zero else -one
                else:
                    xs[i] = x[b, i]
            for i in range(N):
                acc = np.float32(0.0)
                for k in range(indptr[i], indptr[i + 1]):
                    acc += data[k] * xs[indices[k]]
                jx[i] = acc
            for i in range(N):
                y[b, i] = y[b, i] + dtf * (ka * x[b, i] + c0f * jx[i])
            for i in range(N):
                xv = x[b, i] + dta * y[b, i]
                if abs(xv) > one:
                    x[b, i] = one if xv > zero else -one
                    y[b, i] = zero
                else:
                    x[b, i] = xv


def sb(g, B, cfg, rng):
    """methods.sb_traj: x, y ~ U(-0.1, 0.1) float32; c0 = xi * 0.5 / (sigma_J sqrt(N)); a(t) = t/S."""
    S = cfg['S']; dt = cfg['dt']
    c0 = cfg.get('xi', 1.0) * 0.5 / (sigma_J(g) * math.sqrt(g.N))
    x = (rng.random((B, g.N), dtype=np.float32) - 0.5) * 0.2
    y = (rng.random((B, g.N), dtype=np.float32) - 0.5) * 0.2
    _sb_kernel(g.indptr, g.indices, g.data, x, y, S, dt, c0, cfg['family'] == 'dSB')
    return np.where(x >= 0, 1.0, -1.0).astype(np.float32)


# ------------------------------------------------------------------ adiabatic SB (Goto et al. 2019)
@nb.njit(cache=True)
def _asb_kernel(indptr, indices, data, x, y, S, M, delta, kick):
    B, N = x.shape
    jx = np.zeros(N, np.float32)
    one = np.float32(1.0)
    finite = np.ones(B, np.bool_)
    for b in range(B):
        for t in range(1, S + 1):
            p = np.float32(t / S)
            for _ in range(M):
                for i in range(N):
                    x[b, i] = x[b, i] + y[b, i] * delta
                for i in range(N):
                    xv = x[b, i]
                    y[b, i] = y[b, i] - ((xv * xv + one - p) * xv) * delta
            for i in range(N):
                acc = np.float32(0.0)
                for k in range(indptr[i], indptr[i + 1]):
                    acc += data[k] * x[b, indices[k]]
                jx[i] = acc
            for i in range(N):
                y[b, i] = y[b, i] + jx[i] * kick
            ok = True
            for i in range(N):
                if not (np.isfinite(x[b, i]) and np.isfinite(y[b, i])):
                    ok = False
            if not ok:
                finite[b] = False
                for i in range(N):
                    if not np.isfinite(x[b, i]):
                        x[b, i] = 0.0
                    if not np.isfinite(y[b, i]):
                        y[b, i] = 0.0
    return finite


def asb(g, B, cfg, rng):
    """methods.asb_traj with the sigma_J normalisation; returns spins and a per-run finiteness flag."""
    S = cfg['S']; dt = cfg['dt']; M = cfg.get('M', 5)
    xi = np.float32(cfg['xi'] * 0.7 / (sigma_J(g) * math.sqrt(g.N)))
    s0 = rng.choice(np.array([-1.0, 1.0], dtype=np.float32), size=(B, g.N))
    x = np.zeros((B, g.N), np.float32); y = (s0 * 0.1 * rng.random((B, g.N), dtype=np.float32)).astype(np.float32)
    delta = np.float32(dt / M)
    finite = _asb_kernel(g.indptr, g.indices, g.data, x, y, S, M, delta, np.float32(xi * np.float32(dt)))
    return np.where(x >= 0, 1.0, -1.0).astype(np.float32), finite


# ------------------------------------------------------------------ ReAIM ASA (pilot.Batch with sparse field update)
@nb.njit(cache=True)
def _scatter(h, delta, indptr, indices, data):
    B, N = delta.shape
    for b in range(B):
        for i in range(N):
            d = delta[b, i]
            if d != 0.0:
                for k in range(indptr[i], indptr[i + 1]):
                    h[b, indices[k]] += d * data[k]


class Batch:
    """research/reaim_reproduction_20261003/pilot.py::Batch, identical except the exact sparse field update."""

    def __init__(self, g, x, fifo, F=np.max, h=None):
        self.g, self.x, self.fifo, self.F = g, x, fifo, F
        self.h = g.field(x) if h is None else h

    def energy(self):
        return -0.5 * np.einsum("bi,bi->b", self.x, self.h)

    def step(self, k, T, rng):
        x, h = self.x, self.h
        d = 2.0 * x * h
        nN = (d < 0).sum(1)
        q = self.F(self.fifo, axis=1)
        pos = np.where(d > 0, d, np.inf).min(1)
        dmin, dmax = d.min(1), d.max(1)
        th_pm = np.where(np.isfinite(pos), pos, dmax)
        th_mn = dmin + rng.random(len(x)) * T * (dmax - dmin)
        th = np.where(nN > q, th_pm, th_mn)
        C = d <= th[:, None]
        keys = np.where(C, rng.random(C.shape), -1.0)
        kk = np.broadcast_to(np.asarray(k), (len(x),))
        kmax = int(kk.max())
        idx = np.argpartition(-keys, kmax - 1, axis=1)[:, :kmax]
        sel = np.take_along_axis(keys, idx, 1) >= 0
        sel &= np.arange(kmax)[None, :] < kk[:, None]
        rows = np.nonzero(sel)
        flat = idx[rows]
        delta = np.zeros_like(x)
        delta[rows[0], flat] = -2.0 * x[rows[0], flat]
        x += delta
        _scatter(h, delta, self.g.indptr, self.g.indices, self.g.data)
        self.fifo = np.concatenate([self.fifo[:, 1:], nN[:, None]], axis=1)


def reaim(g, B, cfg, rng):
    """methods.reaim_traj (final state of the search trajectory). k values are absolute flip counts."""
    S = cfg['S']; kset = cfg['kset']; t0, t1 = 1.0, cfg['T1']; it_trial, it_run = 32, 96
    N = g.N; R = len(kset)
    alpha = (t1 / t0) ** (1.0 / (S - 1))
    x = rng.choice(np.array([-1.0, 1.0], dtype=np.float32), size=(B, N))
    cur = Batch(g, x, np.full((B, 20), N, dtype=np.int64), np.max)
    t, T = 0, t0
    while t < S:
        n = min(it_trial, S - t)
        rep = Batch(g, np.repeat(cur.x, R, 0), np.repeat(cur.fifo, R, 0), np.max, h=np.repeat(cur.h, R, 0))
        kvec = np.tile(np.array(kset), B)
        Elast = None
        for u in range(n):
            rep.step(kvec, T, rng); T *= alpha
        Elast = rep.energy().reshape(B, R)
        pick = Elast.argmin(1)
        t += n
        sel = np.arange(B) * R + pick
        cur = Batch(g, rep.x[sel].copy(), rep.fifo[sel].copy(), np.max, h=rep.h[sel].copy())
        kbest = np.array(kset)[pick]
        n = min(it_run, S - t)
        for u in range(n):
            cur.step(kbest, T, rng); T *= alpha
        t += n
    return cur.x
