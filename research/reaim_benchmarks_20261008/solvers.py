"""The ten methods on a problems.Problem (CSR couplings + uniform coupling gamma + bias b).

Ports of research/gset_20261007/engine.py and others.py (not edited) with two additions, applied identically to every
method (pre-declared in PROTOCOL.md):
  * bias: the local field is h_i = sum_j J_ij s_j + b_i. Incremental methods (the SCA engine family, SA, ReAIM) put b
    into the fields at initialization and update them per flip, so the bias persists exactly (the hardware scheme of
    MULTIBIT_SPEC.md); SB methods add b_i to every J.x product. Both equal a frozen auxiliary spin s_0 = +1 with
    J_i0 = b_i that never updates and is not counted anywhere (n_lin, statistics).
  * uniform coupling gamma (GPP balance penalty): field += gamma * (M - s_i), M = sum_j s_j tracked exactly.
With gamma = 0 and b = 0 every kernel performs exactly the operations of the original (verify.py checks bit-identity
against engine.run, others.sa, others.reaim and others.sb on G-set instances).
"""
import math

import numpy as np
import numba as nb

ENGINE_FAMILIES = ('plain', 'tec', 'tecT', 'onsager', 'apc')


def ramp_factor(t, S, ramp):            # = research/ablation_20261005/abl.py::ramp_factor
    return 1.0 if (not ramp or t < 0.7 * S) else (S - 1 - t) / max(1.0, S - 1 - 0.7 * S)


def sched(t0, S, tfin):                 # = engine.sched
    return t0 * (tfin / t0) ** (np.arange(S) / max(1, S - 1))


@nb.njit(cache=True)
def field_csr_b(s, indptr, indices, data, b):
    """h = s @ Js + b (float32, exact for integer data): engine.field_csr plus the bias."""
    B, N = s.shape
    h = np.zeros((B, N), np.float32)
    for bb in range(B):
        for i in range(N):
            acc = np.float32(0.0)
            for k in range(indptr[i], indptr[i + 1]):
                acc += data[k] * s[bb, indices[k]]
            h[bb, i] = acc + b[i]
    return h


# ===================================================================== SCA engine family (engine._step + gamma)
@nb.njit(cache=True)
def _step(mode, count, apc, s, h, Mv, gamma, sp, q, u, Tt, c32, c64, cprev, cnt, flips,
          indptr, indices, data, q_reset, r_q, q_lim, flipbuf):
    B, N = s.shape
    d4f = np.float32(4.0 * Tt); lo32 = np.float32(-2.0 * Tt); hi32 = np.float32(2.0 * Tt)
    d464 = 4.0 * Tt; lo64 = -2.0 * Tt; hi64 = 2.0 * Tt
    half = np.float32(0.5); z0 = np.float32(0.0); o1 = np.float32(1.0)
    m2 = np.float32(-2.0)
    g32 = np.float32(gamma); usegam = gamma != 0.0
    for b in range(B):
        nlin = 0; nfl = 0
        cb = c64 * cprev[b]
        Mb = Mv[b]
        for i in range(N):
            hb = h[b, i]
            if usegam:
                hb = hb + g32 * (Mb - s[b, i])
            if mode == 3:
                fld = np.float64(hb) - cb * np.float64(sp[b, i])
                z = np.float64(s[b, i]) * fld + np.float64(q[b, i])
                x = z / d464 + 0.5
                if x < 0.0:
                    x = 0.0
                elif x > 1.0:
                    x = 1.0
                f = x < np.float64(u[b, i])
                if count and z > lo64 and z < hi64:
                    nlin += 1
            else:
                if mode == 0:
                    fld32 = hb
                elif mode == 1:
                    fld32 = hb + c32 * sp[b, i]
                else:
                    fld32 = hb - c32 * sp[b, i]
                z32 = s[b, i] * fld32 + q[b, i]
                x32 = z32 / d4f + half
                if x32 < z0:
                    x32 = z0
                elif x32 > o1:
                    x32 = o1
                f = x32 < u[b, i]
                if count and z32 > lo32 and z32 < hi32:
                    nlin += 1
            if apc:
                if f:
                    q[b, i] = q_reset
                else:
                    qq = q[b, i] * r_q
                    q[b, i] = qq if qq > q_lim else q_lim
            if f:
                flipbuf[nfl] = i; nfl += 1
        cnt[b] = nlin
        flips[b] += nfl
        for i in range(N):
            sp[b, i] = s[b, i]
        for k in range(nfl):
            i = flipbuf[k]
            d = m2 * s[b, i]
            s[b, i] = s[b, i] + d
            Mv[b] += d
            for e in range(indptr[i], indptr[i + 1]):
                h[b, indices[e]] += d * data[e]


def engine_run(P, s0, cfg, rng, tfin, trace_nlin=False):
    """engine.run for a Problem. Families plain, tec, tecT, onsager, apc. Returns spins (B, N) float32, flips (B,)."""
    fam = cfg['family']; S = cfg['S']; T = sched(cfg['T0'], S, tfin)
    s = s0.astype(np.float32).copy(); h = field_csr_b(s, P.indptr, P.indices, P.data, P.b)
    B, N = s.shape
    Mv = s.sum(1).astype(np.float32)
    flips = np.zeros(B); sp = np.zeros_like(s); cprev = np.zeros(B); cnt = np.zeros(B, np.int64)
    q = np.full((B, N), float(cfg.get('q', 0.0)), np.float32)
    apc = fam == 'apc'
    if apc:
        q[:] = cfg['q_reset']
    q_reset = np.float32(cfg.get('q_reset', 0.0)); r_q = np.float32(cfg.get('r_q', 0.0)); q_lim = np.float32(cfg.get('q_lim', 0.0))
    flipbuf = np.zeros(N, np.int64)
    T_prev = None
    nl = np.zeros((S, B), np.int64) if trace_nlin else None
    for t, Tt in enumerate(T):
        mode, c32, c64 = 0, np.float32(0.0), 0.0
        if t > 0:
            if fam == 'onsager':
                mode = 3; c64 = cfg['lam'] * ramp_factor(t, S, cfg['ramp'])
            elif fam == 'tecT':
                mode = 2; c32 = np.float32(cfg['kappa'] * T_prev * ramp_factor(t, S, cfg['ramp']))
            elif fam == 'tec':
                mode = 1; c32 = np.float32(cfg['jv'])
        u = rng.random(s.shape, dtype=np.float32)
        _step(mode, fam == 'onsager' or trace_nlin, apc, s, h, Mv, P.gamma, sp, q, u, float(Tt), c32, float(c64), cprev,
              cnt, flips, P.indptr, P.indices, P.data, q_reset, r_q, q_lim, flipbuf)
        if trace_nlin:
            nl[t] = cnt
        if fam == 'onsager':
            cprev = cnt / (2 * Tt)
        T_prev = Tt
    if trace_nlin:
        return s, flips, nl
    return s, flips


# ===================================================================== SA (others._sa_kernel_csr + bias + gamma)
@nb.njit(cache=True)
def _sa_kernel(indptr, indices, data, bias, gamma, s, Tsched, seeds):
    B, N = s.shape
    S = Tsched.shape[0]
    usegam = gamma != 0.0
    for b in range(B):
        np.random.seed(seeds[b])
        sb = s[b].copy()
        h = np.zeros(N)
        M = 0.0
        for i in range(N):
            acc = 0.0
            for k in range(indptr[i], indptr[i + 1]):
                acc += data[k] * sb[indices[k]]
            h[i] = acc + bias[i]
            M += sb[i]
        for t in range(S):
            Tt = Tsched[t]
            for i in range(N):
                hi = h[i]
                if usegam:
                    hi = hi + gamma * (M - sb[i])
                dE = 2.0 * sb[i] * hi
                if dE <= 0.0 or np.random.random() < math.exp(-dE / Tt):
                    sb[i] = -sb[i]
                    d = 2.0 * sb[i]
                    M += d
                    for k in range(indptr[i], indptr[i + 1]):
                        h[indices[k]] += d * data[k]
        s[b] = sb


def sa(P, B, cfg, rng):
    s0 = rng.choice(np.array([-1.0, 1.0], dtype=np.float32), size=(B, P.N))
    seeds = rng.integers(0, 2 ** 31 - 1, size=B)
    S = cfg['S']
    Tsched = cfg['T0'] * (cfg['T1'] / cfg['T0']) ** (np.arange(S) / max(1, S - 1))
    s = s0.astype(np.float64).copy()
    _sa_kernel(P.indptr, P.indices, P.data.astype(np.float64), P.b.astype(np.float64), float(P.gamma), s,
               Tsched.astype(np.float64), np.asarray(seeds, np.int64))
    return s.astype(np.float32)


# ===================================================================== bSB / dSB (others._sb_kernel + bias + gamma)
@nb.njit(cache=True)
def _sb_kernel(indptr, indices, data, bias, gamma, x, y, S, dt, c0, disc):
    B, N = x.shape
    a0 = 1.0
    jx = np.zeros(N, np.float32); xs = np.zeros(N, np.float32)
    dtf = np.float32(dt); c0f = np.float32(c0); dta = np.float32(dt * a0)
    one = np.float32(1.0); zero = np.float32(0.0)
    g32 = np.float32(gamma); usegam = gamma != 0.0
    for b in range(B):
        for t in range(S):
            at = a0 * t / S
            ka = np.float32(-(a0 - at))
            for i in range(N):
                if disc:
                    xs[i] = one if x[b, i] >= zero else -one
                else:
                    xs[i] = x[b, i]
            X = np.float32(0.0)
            if usegam:
                for i in range(N):
                    X += xs[i]
            for i in range(N):
                acc = np.float32(0.0)
                for k in range(indptr[i], indptr[i + 1]):
                    acc += data[k] * xs[indices[k]]
                if usegam:
                    acc += g32 * (X - xs[i])
                jx[i] = acc + bias[i]
            for i in range(N):
                y[b, i] = y[b, i] + dtf * (ka * x[b, i] + c0f * jx[i])
            for i in range(N):
                xv = x[b, i] + dta * y[b, i]
                if abs(xv) > one:
                    x[b, i] = one if xv > zero else -one
                    y[b, i] = zero
                else:
                    x[b, i] = xv


def sb(P, B, cfg, rng):
    """others.sb: x, y ~ U(-0.1, 0.1) float32; c0 = xi * 0.5 / (sigma_J sqrt(N)); a(t) = t/S."""
    S = cfg['S']; dt = cfg['dt']
    c0 = cfg.get('xi', 1.0) * 0.5 / (P.stats()['sigmaJ_goto'] * math.sqrt(P.N))
    x = (rng.random((B, P.N), dtype=np.float32) - 0.5) * 0.2
    y = (rng.random((B, P.N), dtype=np.float32) - 0.5) * 0.2
    _sb_kernel(P.indptr, P.indices, P.data, P.b, float(P.gamma), x, y, S, dt, c0, cfg['family'] == 'dSB')
    return np.where(x >= 0, 1.0, -1.0).astype(np.float32)


# ===================================================================== aSB (others._asb_kernel + bias + gamma)
@nb.njit(cache=True)
def _asb_kernel(indptr, indices, data, bias, gamma, x, y, S, M, delta, kick):
    B, N = x.shape
    jx = np.zeros(N, np.float32)
    one = np.float32(1.0)
    finite = np.ones(B, np.bool_)
    g32 = np.float32(gamma); usegam = gamma != 0.0
    for b in range(B):
        for t in range(1, S + 1):
            p = np.float32(t / S)
            for _ in range(M):
                for i in range(N):
                    x[b, i] = x[b, i] + y[b, i] * delta
                for i in range(N):
                    xv = x[b, i]
                    y[b, i] = y[b, i] - ((xv * xv + one - p) * xv) * delta
            X = np.float32(0.0)
            if usegam:
                for i in range(N):
                    X += x[b, i]
            for i in range(N):
                acc = np.float32(0.0)
                for k in range(indptr[i], indptr[i + 1]):
                    acc += data[k] * x[b, indices[k]]
                if usegam:
                    acc += g32 * (X - x[b, i])
                jx[i] = acc + bias[i]
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


def asb(P, B, cfg, rng):
    S = cfg['S']; dt = cfg['dt']; M = cfg.get('M', 5)
    xi = np.float32(cfg['xi'] * 0.7 / (P.stats()['sigmaJ_goto'] * math.sqrt(P.N)))
    s0 = rng.choice(np.array([-1.0, 1.0], dtype=np.float32), size=(B, P.N))
    x = np.zeros((B, P.N), np.float32); y = (s0 * 0.1 * rng.random((B, P.N), dtype=np.float32)).astype(np.float32)
    delta = np.float32(dt / M)
    finite = _asb_kernel(P.indptr, P.indices, P.data, P.b, float(P.gamma), x, y, S, M, delta,
                         np.float32(xi * np.float32(dt)))
    return np.where(x >= 0, 1.0, -1.0).astype(np.float32), finite


# ===================================================================== ReAIM ASA (others.Batch + bias + gamma)
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
    """others.Batch with h = Js x + b (incremental) and the full field h + gamma (M - x)."""

    def __init__(self, P, x, fifo, F=np.max, h=None, M=None):
        self.P, self.x, self.fifo, self.F = P, x, fifo, F
        self.h = field_csr_b(np.ascontiguousarray(x, dtype=np.float32), P.indptr, P.indices, P.data, P.b) if h is None else h
        self.M = x.sum(1).astype(np.float32) if M is None else M

    def full(self):
        if self.P.gamma == 0.0:
            return self.h
        return self.h + np.float32(self.P.gamma) * (self.M[:, None] - self.x)

    def energy(self):
        hf = self.full()
        if not np.any(self.P.b):
            return -0.5 * np.einsum("bi,bi->b", self.x, hf)
        return -0.5 * np.einsum("bi,bi->b", self.x, hf + self.P.b[None, :])

    def step(self, k, T, rng):
        x = self.x
        h = self.full()
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
        _scatter(self.h, delta, self.P.indptr, self.P.indices, self.P.data)
        if self.P.gamma != 0.0:
            self.M += delta.sum(1)
        self.fifo = np.concatenate([self.fifo[:, 1:], nN[:, None]], axis=1)


def reaim(P, B, cfg, rng):
    """others.reaim (final state of the search trajectory). k values are absolute flip counts."""
    S = cfg['S']; kset = cfg['kset']; t0, t1 = 1.0, cfg['T1']; it_trial, it_run = 32, 96
    N = P.N; R = len(kset)
    alpha = (t1 / t0) ** (1.0 / (S - 1))
    x = rng.choice(np.array([-1.0, 1.0], dtype=np.float32), size=(B, N))
    cur = Batch(P, x, np.full((B, 20), N, dtype=np.int64), np.max)
    t, T = 0, t0
    while t < S:
        n = min(it_trial, S - t)
        rep = Batch(P, np.repeat(cur.x, R, 0), np.repeat(cur.fifo, R, 0), np.max, h=np.repeat(cur.h, R, 0),
                    M=np.repeat(cur.M, R, 0))
        kvec = np.tile(np.array(kset), B)
        for u in range(n):
            rep.step(kvec, T, rng); T *= alpha
        Elast = rep.energy().reshape(B, R)
        pick = Elast.argmin(1)
        t += n
        sel = np.arange(B) * R + pick
        cur = Batch(P, rep.x[sel].copy(), rep.fifo[sel].copy(), np.max, h=rep.h[sel].copy(), M=rep.M[sel].copy())
        kbest = np.array(kset)[pick]
        n = min(it_run, S - t)
        for u in range(n):
            cur.step(kbest, T, rng); T *= alpha
        t += n
    return cur.x


# ===================================================================== dispatcher
def run_method(P, cfg, B, seed):
    """Same RNG use as research/gset_20261007/run_gset.py::run_method. Returns (spins, flips or None, finite or None)."""
    rng = np.random.default_rng(seed)
    fam = cfg['family']; flips = None; finite = None
    if fam in ENGINE_FAMILIES:
        s0 = rng.choice(np.array([-1.0, 1.0], dtype=np.float32), size=(B, P.N))
        s, flips = engine_run(P, s0, cfg, rng, cfg['tfin'])
    elif fam == 'SA':
        s = sa(P, B, cfg, rng)
    elif fam in ('bSB', 'dSB'):
        s = sb(P, B, cfg, rng)
    elif fam == 'aSB':
        s, finite = asb(P, B, cfg, rng)
    elif fam == 'ReAIM':
        s = reaim(P, B, cfg, rng)
    else:
        raise ValueError(fam)
    return s, flips, finite
