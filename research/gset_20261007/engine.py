"""Sparse, bit-faithful version of the engine model research/ablation_20261005/abl.py::run (not edited) for G-set.

abl.run(J, s0, cfg, rng) is the float model of the Anechoic SCA engine used by every ablation and hardware schedule. Here:
  * the schedule is T0 -> T_fin geometric (abl.run hard-wires T_fin = 5, the K2000 value, through sweeps_abc.T_FIN);
  * the field update h += d @ J is a scatter over the flipped spins' CSR rows. With integer couplings and d in {0, +-2}
    every partial sum is an exact float32 integer, so h is bit-identical to the dense BLAS product;
  * the per-spin decision is one fused numba loop that reproduces numpy 1.26's dtype rules of abl.run exactly:
      plain/tec/tecT/apc: float32 field, z = s*field + q, x = z / f32(4T) + 0.5, clip, flip if x < u  (u ~ rng float32)
      onsager (t >= 1):   float64 field h - (lam*ramp*c_prev)*s_prev, float64 z and x, compared with float64(u)
      n_lin (onsager):    (z > -2T) & (z < 2T), float32 comparison at t = 0 (z float32), float64 afterwards
  * the random stream is identical: one rng.random((B, N), dtype=float32) per step, as in abl.run.
verify_engine.py checks bit-equality (final spins and flip counts) against abl.run itself with abl.m.T_FIN patched.
"""
import math
import os
import sys
from pathlib import Path

os.environ.setdefault('NUMBA_CACHE_DIR', str(Path(__file__).resolve().parent / '.numba_cache'))
sys.dont_write_bytecode = True     # never write .pyc files into the (hashed) research folders we import from

import numpy as np  # noqa: E402
import numba as nb  # noqa: E402

R = Path(str(__import__('pathlib').Path(__file__).resolve().parent / '..'))
for _d in ('theory_ideas_20261003', 'reaim_reproduction_20261003', 'ablation_20261005', 'algorithm_compare_20261007'):
    if str(R / _d) not in sys.path:
        sys.path.insert(0, str(R / _d))
import abl  # noqa: E402
abl.J_RESULTS = R / 'theory_ideas_20261003' / 'J_results.json'

ROOT = Path(__file__).resolve().parent
SIGMA_K2000 = math.sqrt(1999.0)       # RMS local field of K2000 under random spins
DN_K2000 = 1999.0 / 2000.0            # mean squared row norm / N for K2000


# ------------------------------------------------------------------ instances
class Graph:
    def __init__(self, name, n=None, edges=None):
        """Load data/<name> (G-set format, 1-based), or build from edges = (i, j, w) arrays (0-based) for synthetic graphs."""
        if edges is None:
            p = ROOT / 'data' / name
            lines = p.read_text().split('\n')
            n, m = (int(x) for x in lines[0].split())
            a = np.array([ln.split() for ln in lines[1:] if ln.strip()], dtype=np.int64)
            assert len(a) == m
            i, j, w = a[:, 0] - 1, a[:, 1] - 1, a[:, 2]
        else:
            i, j, w = (np.asarray(x, np.int64) for x in edges); m = len(i)
        self.name, self.N, self.m = name, n, m
        self.ei, self.ej, self.w = i, j, w
        self.W = int(w.sum())
        # Ising couplings J = -w (symmetric, zero diagonal), CSR rows
        rows = np.concatenate([i, j]); cols = np.concatenate([j, i]); vals = -np.concatenate([w, w]).astype(np.float32)
        order = np.lexsort((cols, rows))
        rows, cols, vals = rows[order], cols[order], vals[order]
        self.indptr = np.zeros(n + 1, np.int64); np.add.at(self.indptr, rows + 1, 1); self.indptr = np.cumsum(self.indptr)
        self.indices = cols.astype(np.int64); self.data = vals
        self.dbar = 2.0 * m / n                       # (1/N) sum_ij J_ij^2
        self.sigma = math.sqrt(self.dbar)             # RMS field under random spins
        self.alpha = self.sigma / SIGMA_K2000         # field-scale factor relative to K2000
        self.dn = self.dbar / n                       # echo-coefficient factor (mean J^2 * d / N)
        self.lam_factor = self.dn / DN_K2000          # lambda_eff = lambda_K2000 * lam_factor

    def dense(self):
        J = np.zeros((self.N, self.N), np.float32)
        J[self.ei, self.ej] = -self.w; J[self.ej, self.ei] = -self.w
        return J

    def field(self, s):
        """h = s @ J for a batch of spins (B, N), exact (integer couplings), float32."""
        return field_csr(np.ascontiguousarray(s, dtype=np.float32), self.indptr, self.indices, self.data)

    def cut(self, s):
        """Exact integer cut of each row of s: sum over edges of w_ij [s_i != s_j]."""
        si = s[:, self.ei]; sj = s[:, self.ej]
        return ((si != sj).astype(np.int64) * self.w[None, :]).sum(1)


@nb.njit(cache=True)
def field_csr(s, indptr, indices, data):
    B, N = s.shape
    h = np.zeros((B, N), np.float32)
    for b in range(B):
        for i in range(N):
            acc = np.float32(0.0)
            for k in range(indptr[i], indptr[i + 1]):
                acc += data[k] * s[b, indices[k]]
            h[b, i] = acc
    return h


# ------------------------------------------------------------------ schedule (same formula as sweeps_abc.sched)
def sched(t0, S, tfin):
    return t0 * (tfin / t0) ** (np.arange(S) / max(1, S - 1))


# ------------------------------------------------------------------ fused step
@nb.njit(cache=True)
def _step(mode, count, apc, s, h, sp, q, u, Tt, c32, c64, cprev, cnt, flips,
          indptr, indices, data, q_reset, r_q, q_lim, flipbuf):
    """mode 0: field = h; 1: h + c32*sp; 2: h - c32*sp (all float32); 3: float64 field h - (c64*cprev[b])*sp."""
    B, N = s.shape
    d4f = np.float32(4.0 * Tt); lo32 = np.float32(-2.0 * Tt); hi32 = np.float32(2.0 * Tt)
    d464 = 4.0 * Tt; lo64 = -2.0 * Tt; hi64 = 2.0 * Tt
    half = np.float32(0.5); z0 = np.float32(0.0); o1 = np.float32(1.0)
    m2 = np.float32(-2.0)
    for b in range(B):
        nlin = 0; nfl = 0
        cb = c64 * cprev[b]
        for i in range(N):
            if mode == 3:
                fld = np.float64(h[b, i]) - cb * np.float64(sp[b, i])
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
                    fld32 = h[b, i]
                elif mode == 1:
                    fld32 = h[b, i] + c32 * sp[b, i]
                else:
                    fld32 = h[b, i] - c32 * sp[b, i]
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
            for e in range(indptr[i], indptr[i + 1]):
                h[b, indices[e]] += d * data[e]


def run(g, s0, cfg, rng, tfin):
    """abl.run for graph g (sparse) with final temperature tfin. Families: plain, tec, tecT, onsager, apc.
    Returns final spins (B, N) float32 and flips per run (float64), exactly as abl.run."""
    fam = cfg['family']; S = cfg['S']; T = sched(cfg['T0'], S, tfin)
    s = s0.astype(np.float32).copy(); h = g.field(s)
    B, N = s.shape
    flips = np.zeros(B); sp = np.zeros_like(s); cprev = np.zeros(B); cnt = np.zeros(B, np.int64)
    q = np.full((B, N), float(cfg.get('q', 0.0)), np.float32)
    apc = fam == 'apc'
    if apc:
        q[:] = cfg['q_reset']
    q_reset = np.float32(cfg.get('q_reset', 0.0)); r_q = np.float32(cfg.get('r_q', 0.0)); q_lim = np.float32(cfg.get('q_lim', 0.0))
    flipbuf = np.zeros(N, np.int64)
    T_prev = None
    for t, Tt in enumerate(T):
        mode, c32, c64 = 0, np.float32(0.0), 0.0
        if t > 0:
            if fam == 'onsager':
                mode = 3; c64 = cfg['lam'] * abl.ramp_factor(t, S, cfg['ramp'])
            elif fam == 'tecT':
                mode = 2; c32 = np.float32(cfg['kappa'] * T_prev * abl.ramp_factor(t, S, cfg['ramp']))
            elif fam == 'tec':
                mode = 1; c32 = np.float32(cfg['jv'])
        u = rng.random(s.shape, dtype=np.float32)
        _step(mode, fam == 'onsager', apc, s, h, sp, q, u, float(Tt), c32, float(c64), cprev, cnt, flips,
              g.indptr, g.indices, g.data, q_reset, r_q, q_lim, flipbuf)
        if fam == 'onsager':
            cprev = cnt / (2 * Tt)
        T_prev = Tt
    return s, flips


def run_traced(g, s0, cfg, rng, tfin):
    """run() that also returns the mean (over runs) of the applied temporal coefficient per step (diagnostics only)."""
    fam = cfg['family']; S = cfg['S']; T = sched(cfg['T0'], S, tfin)
    s = s0.astype(np.float32).copy(); h = g.field(s)
    B, N = s.shape
    flips = np.zeros(B); sp = np.zeros_like(s); cprev = np.zeros(B); cnt = np.zeros(B, np.int64)
    q = np.full((B, N), float(cfg.get('q', 0.0)), np.float32)
    flipbuf = np.zeros(N, np.int64); T_prev = None; coef = np.zeros(S); band = np.zeros(S)
    for t, Tt in enumerate(T):
        mode, c32, c64 = 0, np.float32(0.0), 0.0
        if t > 0:
            if fam == 'onsager':
                mode = 3; c64 = cfg['lam'] * abl.ramp_factor(t, S, cfg['ramp']); coef[t] = float((c64 * cprev).mean())
            elif fam == 'tecT':
                mode = 2; c32 = np.float32(cfg['kappa'] * T_prev * abl.ramp_factor(t, S, cfg['ramp'])); coef[t] = float(c32)
            elif fam == 'tec':
                mode = 1; c32 = np.float32(cfg['jv']); coef[t] = -float(c32)
        u = rng.random(s.shape, dtype=np.float32)
        _step(mode, True, False, s, h, sp, q, u, float(Tt), c32, float(c64), cprev, cnt, flips,
              g.indptr, g.indices, g.data, np.float32(0), np.float32(0), np.float32(0), flipbuf)
        band[t] = float(cnt.mean()) / N
        if fam == 'onsager':
            cprev = cnt / (2 * Tt)
        T_prev = Tt
    return s, flips, coef, band, T


# ------------------------------------------------------------------ v6.4 cost model (12 x v6.4 engines, 250 MHz)
MHZ = 250.0
E_ENGINES = 12


def cycles(flips, S):
    return 0.1424 * flips + 18.09 * S + 886


def t_trial_ms(flips, S):
    return cycles(flips, S) / (MHZ * 1e3)


def r99(p):
    if p >= 1:
        return 1.0
    if p <= 0:
        return math.inf
    return math.log(0.01) / math.log1p(-p)


def expected_max(x, k):
    """E[max of k i.i.d. draws from the empirical distribution of x] (exact order-statistic formula)."""
    xs = np.sort(np.asarray(x, float)); n = len(xs)
    cdf = np.arange(1, n + 1) / n; prev = np.arange(0, n) / n
    return float((xs * (cdf ** k - prev ** k)).sum())
