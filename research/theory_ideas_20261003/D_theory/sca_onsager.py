"""Light-weight SCA (annealed PCA with inertia q) + Onsager-type lag-2 correction.

Scaled-down (N<=400) version of the N=2000 experiment.  T and q are scaled by
sqrt(N/2000) so that T/sqrt(N) and q/T match the N=2000 schedule.

Update (all spins in parallel, R replicas vectorised):
    h_i   = sum_j J_ij s_j(t)
    corr  = lam * c(t) * s_i(t-1)                    (lag-2 'echo' subtraction)
          [+ optional memory kernel: e(t) = c s(t-1) + gamma e(t-1)]
    z_i   = s_i(t) (h_i - corr_i) + q
    P_stay= clip(z/(4T) + 1/2, 0, 1); flip if P_stay < u
    c(t+1)= n_lin(t) / (2 T_t),  n_lin = #{|z_i| < 2T}
"""
import os
os.environ.setdefault("OMP_NUM_THREADS", "1")
os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")
os.environ.setdefault("VECLIB_MAXIMUM_THREADS", "1")
import numpy as np


def make_J(N, seed, kind="pm1"):
    rng = np.random.default_rng(seed)
    if kind == "pm1":
        A = rng.choice([-1.0, 1.0], size=(N, N))
    else:  # gaussian, unit variance
        A = rng.standard_normal((N, N))
    J = np.triu(A, 1)
    J = J + J.T
    return J


def run(J, lam=0.0, R=64, T0=30.0, T1=5.0, S=560, q=4.0, seed=0, n_ref=2000,
        gamma_mode=None, nlin_uncorrected=False, lam_sched=None, record=False,
        b_target=None):
    """Return dict with final/min energies (per replica) and per-step traces.

    T0,T1,q are given in N=2000 units and rescaled by sqrt(N/n_ref).
    gamma_mode: None (lag-2 only) | 'qchi' (geometric kernel, gamma=q/(2T)).
    lam_sched: optional callable(step, T) -> lam.
    b_target: if set, lam(t) = b_target * 4T^2 / n_lin  (constant 'momentum').
    """
    N = J.shape[0]
    sc = np.sqrt(N / n_ref)
    T0s, T1s, qs = T0 * sc, T1 * sc, q * sc
    rng = np.random.default_rng(seed)
    s = rng.choice([-1.0, 1.0], size=(N, R))
    s_prev = s.copy()
    c = np.zeros(R)
    e = np.zeros((N, R))          # memory-kernel echo estimate
    Ts = T0s * (T1s / T0s) ** (np.arange(S) / (S - 1))
    E_min = np.full(R, np.inf)
    tr = {k: np.zeros(S) for k in ("nlin", "flips", "E", "lam_eff", "x", "q1", "q2", "q3", "fb")} if record else None
    hist = [s.copy()] * 4
    nlin_hist = np.zeros((S, R)) if record else None
    for t in range(S):
        T = Ts[t]
        h = J @ s
        if b_target is not None:
            # lam such that lam*x = b_target with x = n_lin/(4T^2) from previous step
            lam_t = np.where(c > 0, b_target * 2 * T / np.maximum(c, 1e-9), 0.0)
        elif lam_sched is not None:
            lam_t = lam_sched(t, T) * np.ones(R)
        else:
            lam_t = lam * np.ones(R)
        if gamma_mode == "qchi":
            e = c * s_prev + (qs / (2 * T)) * e
            corr = lam_t * e
        else:
            corr = lam_t * c * s_prev
        z = s * (h - corr) + qs
        P = np.clip(z / (4 * T) + 0.5, 0.0, 1.0)
        u = rng.random((N, R))
        flip = P < u
        if nlin_uncorrected:
            z0 = s * h + qs
            nlin = (np.abs(z0) < 2 * T).sum(0)
        else:
            nlin = (np.abs(z) < 2 * T).sum(0)
        c = nlin / (2 * T)
        s_prev = s
        s = np.where(flip, -s, s)
        E = -0.5 * np.einsum("ir,ir->r", s, J @ s)
        E_min = np.minimum(E_min, E)
        if record:
            tr["nlin"][t] = nlin.mean()
            tr["flips"][t] = flip.sum(0).mean()
            tr["E"][t] = E.mean()
            tr["lam_eff"][t] = lam_t.mean()
            tr["x"][t] = (nlin / (4 * T * T)).mean()
            nlin_hist[t] = nlin
            hist = [s] + hist[:3]
            tr["q1"][t] = (hist[0] * hist[1]).mean()
            tr["q2"][t] = (hist[0] * hist[2]).mean()
            tr["q3"][t] = (hist[0] * hist[3]).mean()
            prev_fl = hist[1] != hist[2]          # flipped at previous step
            back = prev_fl & (hist[0] == hist[2])  # ...and flipped back now
            tr["fb"][t] = back.sum() / max(prev_fl.sum(), 1)
    out = {"E_final": E, "E_min": E_min, "Ts": Ts, "q": qs}
    if record:
        out["trace"] = tr
        out["nlin_hist"] = nlin_hist
    return out
