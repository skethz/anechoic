"""Fast CPU model of the Anechoic SCA engine families for long schedules (numba, one trial per thread).

Same algorithm as research/ablation_20261005/abl.py::run (clipped synchronous SCA, Eq. 7):
    z_i   = s_i h_i + q_t - corr_t * s_i * sprev_i
    flip  iff clip(z_i / (4 T_t) + 1/2, 0, 1) < u_i,     u_i ~ U[0, 1)
    n_lin = #{ -2T_t < z_i < 2T_t }
    corr_t = a_t * n_lin(t-1) + b_t        (a_0 = b_0 = 0: no s_prev at t = 0)
  plain     : a = b = 0
  onsager   : a_t = lam * ramp(t) / (2 T_{t-1})          (Onsager-online, coefficient n_lin / 2T)
  tecT      : b_t = kappa * T_{t-1} * ramp(t)            (Onsager-kT)
  ramp(t)   = 1 for t < 0.7 S, then linear to 0 at t = S - 1 (abl.ramp_factor)
  T_t       = T0 * (T1 / T0) ** (t / (S - 1))             (abl uses T1 = 5)
Optional finishing phase (family c): L extra steps at constant (fin_T, fin_q) with no correction (a = b = 0), appended
after the anneal; these are extra rows of the engine's per-step tables (fourT, q, kcorr, kconst).

Differences from abl.run, none of which changes the algorithm:
  * integer spins/fields (int8 / int16) and a sparse field update (only flipped rows are added), instead of float32 BLAS;
  * one splitmix64 stream per trial; a uniform is drawn only when 0 < z/4T + 1/2 < 1 (otherwise the outcome is fixed);
  * z is evaluated in float64.
Validation: test_sca_fast.py (exact replay of abl.run with injected uniforms, and a statistical comparison).
"""
import math
import sys
from pathlib import Path

import numpy as np
import numba as nb

R = Path(str(__import__('pathlib').Path(__file__).resolve().parent / '..'))
for _d in ('theory_ideas_20261003', 'reaim_reproduction_20261003', 'ablation_20261005', 'algorithm_compare_20261007'):
    sys.path.insert(0, str(R / _d))
import abl  # noqa: E402
abl.J_RESULTS = R / 'theory_ideas_20261003' / 'J_results.json'
m = abl.m

N_SPINS = 2000
BEST_KNOWN = 33337

GOLD = np.uint64(0x9E3779B97F4A7C15)
M1 = np.uint64(0xBF58476D1CE4E5B9)
M2 = np.uint64(0x94D049BB133111EB)
S30 = np.uint64(30)
S27 = np.uint64(27)
S31 = np.uint64(31)
S40 = np.uint64(40)
INV24 = 1.0 / 16777216.0


def load():
    """J (float32, J = -W, zero diagonal), sumw, and the int8 copy used by the kernel."""
    J, sumw = m.load()
    J8 = J.astype(np.int8)
    assert np.array_equal(J8.astype(np.float32), J)
    return J, float(sumw), J8


# ----------------------------------------------------------------------------------------------- schedule tables
def ramp_vec(S, ramp, start=0.7):
    t = np.arange(S, dtype=np.float64)
    if not ramp:
        return np.ones(S)
    r = (S - 1 - t) / max(1.0, S - 1 - start * S)
    return np.where(t < start * S, 1.0, r)


def tables(cfg):
    """Per-step tables (T, q, a, b) of length S + fin_L for a configuration dict."""
    fam = cfg['family']; S = int(cfg['S'])
    T0 = float(cfg['T0']); T1 = float(cfg.get('T1', 5.0))
    T = T0 * (T1 / T0) ** (np.arange(S) / max(1, S - 1))
    q = np.full(S, float(cfg.get('q', 0.0)))
    a = np.zeros(S); b = np.zeros(S)
    rp = ramp_vec(S, cfg.get('ramp', False), cfg.get('ramp_start', 0.7))
    if fam == 'onsager':
        a[1:] = cfg['lam'] * rp[1:] / (2.0 * T[:-1])
    elif fam == 'tecT':
        b[1:] = cfg['kappa'] * T[:-1] * rp[1:]
    elif fam == 'tec':
        b[1:] = -cfg['jv']
    elif fam != 'plain':
        raise ValueError(fam)
    L = int(cfg.get('fin_L', 0))
    if L > 0:
        T = np.concatenate([T, np.full(L, float(cfg['fin_T']))])
        q = np.concatenate([q, np.full(L, float(cfg['fin_q']))])
        a = np.concatenate([a, np.zeros(L)]); b = np.concatenate([b, np.zeros(L)])
    return T.astype(np.float64), q.astype(np.float64), a.astype(np.float64), b.astype(np.float64)


# ----------------------------------------------------------------------------------------------- kernel
@nb.njit(inline='always')
def _next(x):
    x = x + GOLD
    z = x
    z = (z ^ (z >> S30)) * M1
    z = (z ^ (z >> S27)) * M2
    z = z ^ (z >> S31)
    return x, z


@nb.njit(parallel=True, cache=True, fastmath=False)
def run_steps(J2, s, sp, h, nlin, rng, flips, best2, Ttab, qtab, atab, btab, t0, t1, U, step_flips):
    """Advance every trial b through steps t0 .. t1-1 in place.
    J2: int8 (N, N) = 2*J.  s, sp: int8 (B, N).  h: int16 (B, N) = J s.  nlin, flips, best2: int64 (B,).
    rng: uint64 (B,) splitmix64 states.  best2: running max of sum_i s_i h_i over visited states (cut = (sumw + best2/2)/2).
    U: float32 (S, B, N) injected uniforms for the exactness test, or shape (0, 0, 0) for the internal RNG.
    step_flips: int32 (B, S) per-step flip counts, or shape (0, 0) to skip."""
    B, N = s.shape
    useU = U.shape[0] > 0
    rec = step_flips.shape[0] > 0
    for b in nb.prange(B):
        x = rng[b]
        buf = np.empty(N, np.int32)
        sb = s[b]; spb = sp[b]; hb = h[b]
        nl = nlin[b]; fl = flips[b]; bst = best2[b]
        for t in range(t0, t1):
            T = Ttab[t]; fourT = 4.0 * T; twoT = 2.0 * T
            corr = atab[t] * nl + btab[t]
            q = qtab[t]
            nf = 0; cnt = 0; sh = 0
            for i in range(N):
                si = np.int64(sb[i])
                shi = si * np.int64(hb[i])
                sh += shi
                if sb[i] == spb[i]:
                    z = shi + q - corr
                else:
                    z = shi + q + corr
                if z <= -twoT:
                    if useU:
                        if U[t, b, i] > 0.0:
                            buf[nf] = i; nf += 1
                    else:
                        buf[nf] = i; nf += 1
                elif z < twoT:
                    cnt += 1
                    if useU:
                        u = np.float64(U[t, b, i])
                    else:
                        x, r = _next(x)
                        u = np.float64(r >> S40) * INV24
                    if z / fourT + 0.5 < u:
                        buf[nf] = i; nf += 1
            if sh > bst:
                bst = sh
            for i in range(N):
                spb[i] = sb[i]
            for k in range(nf):
                xx = buf[k]
                sb[xx] = -sb[xx]
            for k in range(nf):
                xx = buf[k]
                row = J2[xx]
                if sb[xx] > 0:
                    for y in range(N):
                        hb[y] += row[y]
                else:
                    for y in range(N):
                        hb[y] -= row[y]
            fl += nf
            nl = cnt
            if rec:
                step_flips[b, t] = nf
        sh = 0
        for i in range(N):
            sh += np.int64(sb[i]) * np.int64(hb[i])
        if sh > bst:
            bst = sh
        rng[b] = x; nlin[b] = nl; flips[b] = fl; best2[b] = bst


@nb.njit(parallel=True, cache=True)
def sh_sum(s, h):
    B, N = s.shape
    out = np.empty(B, np.int64)
    for b in nb.prange(B):
        acc = 0
        for i in range(N):
            acc += np.int64(s[b, i]) * np.int64(h[b, i])
        out[b] = acc
    return out


@nb.njit(cache=True)
def _greedy_one(J2, sb, hb):
    N = sb.shape[0]
    n = 0
    while True:
        best = 0
        bi = -1
        for i in range(N):
            v = np.int64(sb[i]) * np.int64(hb[i])
            if v < best:
                best = v
                bi = i
        if bi < 0:
            return n
        sb[bi] = -sb[bi]
        row = J2[bi]
        if sb[bi] > 0:
            for y in range(N):
                hb[y] += row[y]
        else:
            for y in range(N):
                hb[y] -= row[y]
        n += 1


@nb.njit(parallel=True, cache=True)
def greedy_descent(J2, s, h):
    """Steepest single-flip descent to a 1-flip local optimum, in place (host-side reference for family c).
    Repeatedly flips argmin_i s_i h_i while it is negative (cut gain = -s_i h_i). Returns flips per trial."""
    B = s.shape[0]
    nfl = np.zeros(B, np.int64)
    for b in nb.prange(B):
        nfl[b] = _greedy_one(J2, s[b], h[b])
    return nfl


# ----------------------------------------------------------------------------------------------- state helpers
class State:
    """Batch of B trials (or replicas)."""

    def __init__(self, J, J8, seeds_or_s0, rng_seeds=None):
        if isinstance(seeds_or_s0, np.ndarray) and seeds_or_s0.ndim == 2:
            s0 = seeds_or_s0.astype(np.int8)
        else:
            raise TypeError('pass s0 (B, N) array')
        self.s = np.ascontiguousarray(s0)
        self.sp = self.s.copy()
        hf = self.s.astype(np.float32) @ J
        self.h = np.ascontiguousarray(np.rint(hf).astype(np.int16))
        B = len(self.s)
        self.nlin = np.zeros(B, np.int64)
        self.flips = np.zeros(B, np.int64)
        self.best2 = np.full(B, np.iinfo(np.int64).min // 4, np.int64)
        self.rng = np.asarray(rng_seeds, np.uint64).copy()
        assert self.rng.shape == (B,)


def make_batch(J, B, seed):
    """Independent random initial spins and RNG states for B trials from one integer/SeedSequence seed."""
    ss = np.random.SeedSequence(seed)
    g = np.random.default_rng(ss)
    s0 = np.where(g.random((B, J.shape[0])) < 0.5, -1, 1).astype(np.int8)
    rs = g.integers(0, 2 ** 63 - 1, size=B, dtype=np.int64).astype(np.uint64)
    return s0, rs


_EMPTY_U = np.zeros((0, 0, 0), np.float32)
_EMPTY_SF = np.zeros((0, 0), np.int32)


def advance(J2, st, tabs, t0, t1, U=None, step_flips=None):
    T, q, a, b = tabs
    run_steps(J2, st.s, st.sp, st.h, st.nlin, st.rng, st.flips, st.best2, T, q, a, b, int(t0), int(t1),
              _EMPTY_U if U is None else U, _EMPTY_SF if step_flips is None else step_flips)


def cuts(st, sumw):
    return (sumw + 0.5 * sh_sum(st.s, st.h)) / 2.0


def best_cuts(st, sumw):
    return (sumw + 0.5 * st.best2) / 2.0


def run_cfg(J, J2, sumw, cfg, B, seed, greedy=False, record_steps=False):
    """Run B independent trials of cfg. Returns dict of per-trial arrays."""
    tabs = tables(cfg)
    Stot = len(tabs[0])
    s0, rs = make_batch(J, B, seed)
    st = State(J, None, s0, rs)
    sf = np.zeros((B, Stot), np.int32) if record_steps else None
    advance(J2, st, tabs, 0, Stot, step_flips=sf)
    out = dict(cut=cuts(st, sumw), best=best_cuts(st, sumw), flips=st.flips.copy())
    if sf is not None:
        out['step_flips'] = sf
    if greedy:
        s = st.s.copy(); h = st.h.copy()
        g = greedy_descent(J2, s, h)
        out['cut_greedy'] = (sumw + 0.5 * sh_sum(s, h)) / 2.0
        out['greedy_flips'] = g
    return out


# ----------------------------------------------------------------------------------------------- cost model and statistics
F_V80_MHZ = 250.0
ASIC_GHZ = (1.0, 1.25)
GBSB_TTS_MS = 9.61


def cycles(flips, S_total):
    """v6.4 engine, cycles per trial (measured model): 0.1424 * flips + 18.09 * S + 886."""
    return 0.1424 * np.asarray(flips, dtype=np.float64) + 18.09 * S_total + 886.0


def round_cycles(trial_cycles, E=12, n_boot=4000, seed=12345):
    """Expected max over E concurrent engines of per-trial cycles (bootstrap from the empirical distribution)."""
    c = np.asarray(trial_cycles, np.float64)
    g = np.random.default_rng(seed)
    return float(g.choice(c, size=(n_boot, E), replace=True).max(1).mean())


def r99(p):
    if p >= 1:
        return 1.0
    if p <= 0:
        return math.inf
    return math.log(0.01) / math.log1p(-p)


def tts_ms(t_round_ms, P_round):
    """Eq. 9: TTS99 = t_R * ln(0.01) / ln(1 - P_R); P_R >= 0.99 gives one round."""
    if P_round >= 0.99:
        return t_round_ms
    if P_round <= 0:
        return math.inf
    return t_round_ms * math.log(0.01) / math.log1p(-P_round)


def wilson(k, n, z=1.959963984540054):
    if n == 0:
        return (0.0, 1.0)
    p = k / n; d = 1 + z * z / n
    c = (p + z * z / (2 * n)) / d
    w = z * math.sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / d
    return (max(0.0, c - w), min(1.0, c + w))
