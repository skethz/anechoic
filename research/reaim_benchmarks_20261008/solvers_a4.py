"""Addendum 4 kernels.

  tec_seq  TEC (Du et al. 2026) as a single p-bit network with Glauber dynamics, Eq. (4):
           P_flip,i(t) = 1 / (1 + exp{2 beta sigma_i(t) [h_i^s(t) + J_v sigma_i(t-1)]}), h_i^s = sum_j J_ij sigma_j + b_i,
           sequential single-spin updates in index order within each time step (one sweep = one cycle), the spatial
           field always current (spins relax within the step), the temporal field from the configuration at the end of
           the previous cycle (sigma(t-1); the initial configuration for the first cycle). T geometric T0 -> Tfin over
           the cycles. Optional per-cycle energy trace (validation).
  sca_sig_q  The SCA of Okonogi et al. 2023, Algorithm 1: synchronous update, exact sigmoid flip probability,
           exponential T (Eq. 3) and exponential q(s) = q_init r_q^(s-1) (Eq. 14), output argmin_s H(sigma(s)).
Everything else is dispatched to solvers_a3b.run_method (papers' output rules for APC-SCA and ReAIM; final state
otherwise)."""
import math

import numpy as np
import numba as nb

import solvers as SV
import solvers_a3 as SV3
import solvers_a3b as SV3B


@nb.njit(cache=True)
def _tec_seq_kernel(indptr, indices, data, bias, gamma, s, Tsched, jv, seeds, trace, Etr, flips):
    B, N = s.shape
    S = Tsched.shape[0]
    usegam = gamma != 0.0
    for b in range(B):
        np.random.seed(seeds[b])
        sb = s[b].copy()
        prev = sb.copy()
        h = np.zeros(N)
        M = 0.0
        for i in range(N):
            acc = 0.0
            for k in range(indptr[i], indptr[i + 1]):
                acc += data[k] * sb[indices[k]]
            h[i] = acc + bias[i]
            M += sb[i]
        E = 0.0
        if trace:
            for i in range(N):
                F = h[i]
                if usegam:
                    F += gamma * (M - sb[i])
                E += -0.5 * sb[i] * (F - bias[i]) - bias[i] * sb[i]
        nfl = 0
        for t in range(S):
            beta = 1.0 / Tsched[t]
            for i in range(N):
                hs = h[i]
                if usegam:
                    hs += gamma * (M - sb[i])
                a = 2.0 * beta * sb[i] * (hs + jv * prev[i])
                if a > 700.0:
                    pf = 0.0
                elif a < -700.0:
                    pf = 1.0
                else:
                    pf = 1.0 / (1.0 + math.exp(a))
                if np.random.random() < pf:
                    if trace:
                        E += 2.0 * sb[i] * hs
                    sb[i] = -sb[i]
                    d = 2.0 * sb[i]
                    M += d
                    nfl += 1
                    for k in range(indptr[i], indptr[i + 1]):
                        h[indices[k]] += d * data[k]
            for i in range(N):
                prev[i] = sb[i]
            if trace:
                Etr[b, t] = E
        s[b] = sb
        flips[b] = nfl


def tec_seq(P, B, cfg, rng, trace=False):
    """Random +-1 start (rng.choice), per-run seeds (rng.integers), as solvers.sa. Returns spins, flips, energy trace."""
    s0 = rng.choice(np.array([-1.0, 1.0], dtype=np.float32), size=(B, P.N))
    seeds = rng.integers(0, 2 ** 31 - 1, size=B)
    S = cfg['S']
    T = SV.sched(cfg['T0'], S, cfg['tfin'])
    s = s0.astype(np.float64).copy()
    Etr = np.zeros((B, S) if trace else (1, 1))
    flips = np.zeros(B, np.int64)
    _tec_seq_kernel(P.indptr, P.indices, P.data.astype(np.float64), P.b.astype(np.float64), float(P.gamma), s,
                    T.astype(np.float64), float(cfg['jv']), np.asarray(seeds, np.int64), trace, Etr, flips)
    return s.astype(np.float32), flips.astype(np.float64), (Etr if trace else None)


def sca_sig_q(P, B, cfg, rng):
    """Okonogi 2023 Algorithm 1 (fine-tuned SCA): returns (argmin-energy states, final states, flips)."""
    s0 = rng.choice(np.array([-1.0, 1.0], dtype=np.float32), size=(B, P.N))
    S = cfg['S']; T = SV.sched(cfg['T0'], S, cfg['tfin'])
    qs = cfg['q_init'] * (cfg['q_final'] / cfg['q_init']) ** (np.arange(S) / max(1, S - 1))
    s = s0.astype(np.float64).copy()
    h = SV.field_csr_b(s.astype(np.float32), P.indptr, P.indices, P.data, P.b)
    Mv = s.sum(1)
    flips = np.zeros(B); sp = s.copy(); q = np.zeros((B, P.N))
    flipbuf = np.zeros(P.N, np.int64)
    b64 = P.b.astype(np.float64)
    bestE = SV3B._energy(s, h, Mv, float(P.gamma), b64); best = s.copy()
    for t, Tt in enumerate(T):
        q[:] = qs[t]
        u = rng.random(s.shape)
        SV3._step_sig(False, False, s, h, Mv, float(P.gamma), sp, q, u, float(Tt), 1.0, 0.0, flips,
                      P.indptr, P.indices, P.data, 0.0, 0.0, 0.0, flipbuf)
        E = SV3B._energy(s, h, Mv, float(P.gamma), b64)
        imp = E < bestE
        if imp.any():
            best[imp] = s[imp]; bestE[imp] = E[imp]
    return best.astype(np.float32), s.astype(np.float32), flips


def run_method(P, cfg, B, seed):
    """solvers_a3b.run_method plus 'tec_seq'."""
    if cfg['family'] == 'tec_seq':
        rng = np.random.default_rng(seed)
        s, flips, _ = tec_seq(P, B, cfg, rng)
        return s, flips, None
    return SV3B.run_method(P, cfg, B, seed)
