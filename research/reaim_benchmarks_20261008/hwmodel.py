"""Numba port of the hardware golden reference fpga/v80_sca/src/sca_ref_bias.hpp::run_trial_bias (with sca_ref.hpp's
Rng and make_tables, and the schedule-to-table rules of fpga/v80_sca/src/mb_common.hpp::tables), for the coupling-
precision study. Integer arithmetic throughout: fields h (exact integers), z = s*h*2^16 + q - (s == s_prev ? corr : -corr)
in int64, Q16.16 tables, xoshiro128** per lane (LANES = 256, as built with -DSCA_LANES=256), NP = 2048 slots per step.

Couplings come from a problems.Problem (CSR integer part + uniform gamma + bias): the field h_i = hs_i + gamma (M - s_i)
equals the reference's dense-row field exactly, and the random streams are consumed exactly as in the reference
(2048 draws per step, independent of n). verify_hw.py checks spins, flips, the n_lin trace and sum s*h bit for bit
against the C++ reference on gpu-host.
"""
import math

import numpy as np
import numba as nb

NP = 2048
LANES = 256
ROUNDS = NP // LANES
M32 = np.uint64(0xFFFFFFFF)


def llround(x):
    """std::llround: round half away from zero."""
    return int(math.floor(x + 0.5)) if x >= 0 else -int(math.floor(-x + 0.5))


def geom(t0, t1, S):
    """mb_common.hpp::geom (std::pow in double)."""
    return [t0 * math.pow(t1 / t0, float(t) / max(1, S - 1)) for t in range(S)]


def tables(sch):
    """mb_common.hpp::tables(Sched) with scale = 1: fourT, q (Q16.16), kcorr (Q8.24 of lam_t / 2T_{t-1}), kconst (Q16.16).
    sch: dict(t0, t1, S, q, lam, jv, kappa, ramp)."""
    S = int(sch['S']); T = geom(float(sch['t0']), float(sch['t1']), S)
    lam = float(sch.get('lam', 0.0)); jv = float(sch.get('jv', 0.0)); kappa = float(sch.get('kappa', 0.0))
    ramp = bool(sch.get('ramp', False)); q = float(sch['q'])
    L = [lam] * S
    if ramp:
        for t in range(S):
            if t >= 0.7 * S:
                L[t] = lam * (S - 1 - t) / max(1.0, S - 1 - 0.7 * S)
    fourT = np.array([llround(4.0 * T[t] * 65536.0) for t in range(S)], np.int64)
    qq = np.array([llround(q * 65536.0) for t in range(S)], np.int64)
    kcorr = np.array([0 if t == 0 else llround(L[t] / (2.0 * T[t - 1]) * 16777216.0) for t in range(S)], np.int64)
    kconst = np.array([0 if t == 0 else llround(-jv * 65536.0) for t in range(S)], np.int64)
    if kappa > 0:
        for t in range(1, S):
            rf = (S - 1 - t) / max(1.0, S - 1 - 0.7 * S) if (ramp and t >= 0.7 * S) else 1.0
            kconst[t] = llround(kappa * T[t - 1] * rf * 65536.0)
            kcorr[t] = 0
    return fourT, qq, kcorr, kconst


def table_checks(sch, K=None):
    """Range checks of the golden reference (int32 fourT, q, kconst) and of the FPGA engine
    (mb_common.hpp::check_tables: kcorr fits int32, 0 < 4T < 2^27 in Q16.16, |q| + |corr|max < 2^29 for K <= 4)."""
    fourT, qq, kc, kk = tables(sch)
    i32 = lambda a: bool(np.all(np.abs(a) < 2 ** 31))  # noqa: E731
    corr = ((2048 * np.abs(kc)) >> 8) + np.abs(kk)
    return dict(ref_int32_tables=i32(fourT) and i32(qq) and i32(kk) and i32(kc),
                engine_kcorr_int32=i32(kc), engine_4T_27bit=bool(np.all((fourT > 0) & (fourT < 2 ** 27))),
                engine_q_corr_29bit=bool(K is None or K > 4 or np.all(np.abs(qq) + corr < 2 ** 29)),
                max_T=float(fourT.max() / 4 / 65536.0))


@nb.njit(cache=True)
def _splitmix(x):
    x = x + np.uint64(0x9E3779B97F4A7C15)
    z = x
    z = (z ^ (z >> np.uint64(30))) * np.uint64(0xBF58476D1CE4E5B9)
    z = (z ^ (z >> np.uint64(27))) * np.uint64(0x94D049BB133111EB)
    return x, z ^ (z >> np.uint64(31))


@nb.njit(cache=True)
def _seed_lanes(seed, trial, st):
    for l in range(LANES):
        x = np.uint64(seed) ^ (np.uint64(trial) << np.uint64(20)) ^ (np.uint64(l) << np.uint64(48)) ^ np.uint64(0x5CA0F00D)
        x, a = _splitmix(x)
        x, b = _splitmix(x)
        st[l, 0] = a & M32; st[l, 1] = a >> np.uint64(32); st[l, 2] = b & M32; st[l, 3] = b >> np.uint64(32)
        if (st[l, 0] | st[l, 1] | st[l, 2] | st[l, 3]) == np.uint64(0):
            st[l, 0] = np.uint64(1)


@nb.njit(cache=True)
def _rotl(x, k):
    return ((x << np.uint64(k)) | (x >> np.uint64(32 - k))) & M32


@nb.njit(cache=True)
def _next(st, l):
    s0 = st[l, 0]; s1 = st[l, 1]; s2 = st[l, 2]; s3 = st[l, 3]
    r = (_rotl((s1 * np.uint64(5)) & M32, 7) * np.uint64(9)) & M32
    t = (s1 << np.uint64(9)) & M32
    s2 ^= s0; s3 ^= s1; s1 ^= s2; s0 ^= s3; s2 ^= t; s3 = _rotl(s3, 11)
    st[l, 0] = s0; st[l, 1] = s1; st[l, 2] = s2; st[l, 3] = s3
    return r


@nb.njit(cache=True)
def _trial(n, indptr, indices, data, gamma, bias, fourT, qtab, kcorr, kconst, seed, trial, s_out, nlin_out):
    """One trial; returns (flips, sum_sh, sum_bs). s_out (n,) int64 receives the final spins; nlin_out (S,) the trace."""
    S = fourT.shape[0]
    st = np.zeros((LANES, 4), np.uint64)
    _seed_lanes(seed, trial, st)
    s = np.ones(NP, np.int64); sp = np.ones(NP, np.int64)
    for k in range(ROUNDS):
        for l in range(LANES):
            i = k * LANES + l
            r = _next(st, l)
            if i < n:
                s[i] = 1 if (r >> np.uint64(31)) != np.uint64(0) else -1
    hs = np.zeros(n, np.int64)
    M = 0
    for i in range(n):
        acc = bias[i]
        for e in range(indptr[i], indptr[i + 1]):
            acc += data[e] * s[indices[e]]
        hs[i] = acc
        M += s[i]
    flipped = np.zeros(n, np.int64)
    flips = 0
    nlin_prev = 0
    for t in range(S):
        corr = ((nlin_prev * kcorr[t]) >> 8) + kconst[t]
        f4 = fourT[t]; twoT = f4 >> 1; qt = qtab[t]
        nlin = 0; nf = 0
        for k in range(ROUNDS):
            for l in range(LANES):
                i = k * LANES + l
                u16 = np.int64(_next(st, l) >> np.uint64(16))
                if i >= n:
                    continue
                hi = hs[i] + gamma * (M - s[i])
                c = corr if s[i] == sp[i] else -corr
                z = (s[i] * hi) * 65536 + qt - c
                rr = ((u16 - 32768) * f4) >> 16
                if z > -twoT and z < twoT:
                    nlin += 1
                if z < rr:
                    flipped[nf] = i; nf += 1
        for i in range(n):
            sp[i] = s[i]
        for k in range(nf):
            x = flipped[k]
            s[x] = -s[x]
        for k in range(nf):
            x = flipped[k]
            d = 2 * s[x]
            M += d
            for e in range(indptr[x], indptr[x + 1]):
                hs[indices[e]] += d * data[e]
        flips += nf
        nlin_prev = nlin
        nlin_out[t] = nlin
    sum_sh = 0; sum_bs = 0
    for i in range(n):
        hi = hs[i] + gamma * (M - s[i])
        sum_sh += s[i] * hi
        sum_bs += bias[i] * s[i]
        s_out[i] = s[i]
    return flips, sum_sh, sum_bs


@nb.njit(cache=True)
def _trials(n, indptr, indices, data, gamma, bias, fourT, qtab, kcorr, kconst, seed, trials, S_out, F_out, SH_out, NL_out):
    nlin = np.zeros(fourT.shape[0], np.int64)
    srow = np.zeros(n, np.int64)
    for a in range(trials.shape[0]):
        f, sh, bs = _trial(n, indptr, indices, data, gamma, bias, fourT, qtab, kcorr, kconst, seed, trials[a], srow, nlin)
        S_out[a] = srow
        F_out[a] = f
        SH_out[a, 0] = sh; SH_out[a, 1] = bs
        NL_out[a] = nlin.sum()


def run(P, sch, seed, trials, trace=False):
    """Run trials (array of uint32 trial indices) of the reference arithmetic on Problem P (integer couplings).
    Returns spins (B, n) int8, flips (B,), sum_sh (B,), sum_bs (B,), and the n_lin trace of the last trial if trace."""
    assert P.N <= NP
    fourT, qq, kc, kk = tables(sch)
    data = P.data.astype(np.int64); assert np.all(data.astype(np.float32) == P.data)
    bias = P.b.astype(np.int64); gamma = int(P.gamma); assert gamma == P.gamma
    trials = np.asarray(trials, np.int64)
    B = len(trials)
    So = np.zeros((B, P.N), np.int64); Fo = np.zeros(B, np.int64); SH = np.zeros((B, 2), np.int64); NL = np.zeros(B, np.int64)
    _trials(P.N, P.indptr, P.indices, data, gamma, bias, fourT, qq, kc, kk, np.uint64(seed), trials, So, Fo, SH, NL)
    out = dict(spins=So.astype(np.int8), flips=Fo, sum_sh=SH[:, 0], sum_bs=SH[:, 1], nlin_total=NL)
    if trace:
        nl = np.zeros(int(sch['S']), np.int64); srow = np.zeros(P.N, np.int64)
        _trial(P.N, P.indptr, P.indices, data, gamma, bias, fourT, qq, kc, kk, np.uint64(seed), int(trials[-1]), srow, nl)
        out['nlin_trace_last'] = nl
    return out
