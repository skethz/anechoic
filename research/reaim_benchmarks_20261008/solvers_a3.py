"""Addendum 3 kernels: baselines at their papers' own settings.

  * engine_sig: synchronous (all-spins-at-once) update with the EXACT logistic flip probability, as written in
      - APC-SCA (Okonogi et al. 2023, Algorithm 2): p_i = sigmoid(-(h_i s_i + q_i)/T), per-spin q_i reset to q_reset
        after a flip, else q_i <- max(q_i r_q, q_lim);                                  (beta factor c = 1)
      - TEC (Du et al. 2026, Eq. 4): P_i = 1/(1 + exp{2 beta s_i(t) [h_i(t) + J_v s_i(t-1)]}), no pinning q.  (c = 2)
    Fields h_i = sum_j J_ij s_j + b_i (+ gamma (M - s_i)), exactly as solvers.engine_run. Flip if u < P_flip,
    u ~ U[0, 1) float64 (one rng.random((B, N)) per step). T(t) = T0 (Tfin/T0)^(t/(S-1)) (APC Eq. 3).
  * Everything else is dispatched to solvers_a2.run_method (frozen solvers.py; ReAIM with F and T0).
"""
import math

import numpy as np
import numba as nb

import solvers as SV
import solvers_a2 as SA2


@nb.njit(cache=True)
def _step_sig(apc, tec, s, h, Mv, gamma, sp, q, u, Tt, cb, jv, flips, indptr, indices, data, q_reset, r_q, q_lim, flipbuf):
    B, N = s.shape
    usegam = gamma != 0.0
    for b in range(B):
        nfl = 0
        Mb = Mv[b]
        for i in range(N):
            fld = np.float64(h[b, i])
            if usegam:
                fld += gamma * (Mb - s[b, i])
            if tec:
                fld += jv * sp[b, i]
            z = s[b, i] * fld + q[b, i]
            a = cb * z / Tt
            if a > 700.0:
                pf = 0.0
            elif a < -700.0:
                pf = 1.0
            else:
                pf = 1.0 / (1.0 + math.exp(a))
            f = u[b, i] < pf
            if apc:
                if f:
                    q[b, i] = q_reset
                else:
                    qq = q[b, i] * r_q
                    q[b, i] = qq if qq > q_lim else q_lim
            if f:
                flipbuf[nfl] = i; nfl += 1
        flips[b] += nfl
        for i in range(N):
            sp[b, i] = s[b, i]
        for k in range(nfl):
            i = flipbuf[k]
            d = -2.0 * s[b, i]
            s[b, i] = s[b, i] + d
            Mv[b] += d
            for e in range(indptr[i], indptr[i + 1]):
                h[b, indices[e]] += np.float32(d) * data[e]


def engine_sig(P, s0, cfg, rng):
    """cfg: family 'apc_sig' (q_reset, r_q, q_lim) or 'tec_sig' (jv); T0, tfin, S. Returns spins (B, N), flips (B,)."""
    S = cfg['S']; T = SV.sched(cfg['T0'], S, cfg['tfin'])
    s = s0.astype(np.float64).copy()
    h = SV.field_csr_b(s.astype(np.float32), P.indptr, P.indices, P.data, P.b)
    B, N = s.shape
    Mv = s.sum(1)
    flips = np.zeros(B); sp = s.copy()                  # s(t-1) for t = 0: the initial state (no TEC term difference)
    apc = cfg['family'] == 'apc_sig'; tec = cfg['family'] == 'tec_sig'
    q = np.full((B, N), float(cfg['q_reset']) if apc else float(cfg.get('q', 0.0)))
    flipbuf = np.zeros(N, np.int64)
    for t, Tt in enumerate(T):
        u = rng.random(s.shape)
        jv = float(cfg['jv']) if (tec and t > 0) else 0.0
        _step_sig(apc, tec and t > 0, s, h, Mv, float(P.gamma), sp, q, u, float(Tt), float(cfg['cb']), jv, flips,
                  P.indptr, P.indices, P.data, float(cfg.get('q_reset', 0.0)), float(cfg.get('r_q', 0.0)),
                  float(cfg.get('q_lim', 0.0)), flipbuf)
    return s.astype(np.float32), flips


def run_method(P, cfg, B, seed):
    """solvers_a2.run_method plus the exact-sigmoid families (same RNG creation; random +-1 initial spins)."""
    if cfg['family'] in ('apc_sig', 'tec_sig'):
        rng = np.random.default_rng(seed)
        s0 = rng.choice(np.array([-1.0, 1.0], dtype=np.float32), size=(B, P.N))
        s, flips = engine_sig(P, s0, cfg, rng)
        return s, flips, None
    return SA2.run_method(P, cfg, B, seed)
