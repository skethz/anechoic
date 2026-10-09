"""Dynamical mean-field theory (DMFT) of synchronous stochastic annealers (SCA family) on SK-type couplings.

Rescaled units: A = J/sqrt(N), T~ = T/sqrt(N), q~ = q/sqrt(N), beta~ = beta/sqrt(N).  Raw hardware units are recovered by *sqrt(N).
Real system:  u_i(t) = sum_j A_ij s_j(t);  g_i(t) = u_i(t) + q_t s_i(t) - beta_t s_i(t-1);  P(s_i(t+1)=+1) = (1 + f_t(g))/2
              f_t(g) = clip(g/(2T_t), -1, 1) (clipped, STATICA Eq. 7)  or  tanh(g/(2T_t)) (logistic).
Effective single-site process (N -> inf, symmetric couplings of variance 1/N):
              u(t) = eta(t) + sum_{t'<t} G(t,t') s(t'),   eta Gaussian with E[eta(t)eta(t')] = C(t,t') = E[s(t)s(t')],
              G(t,t') = dE[s(t)]/dtheta(t')  (theta(t') added to g(t')),   C and G self-consistent.
Causality lets one forward pass of M sampled paths solve the equations: at time t, C(t,.) and G(t,.) depend only on s(<=t)
and on the scores at decisions < t. G is estimated with the exact likelihood-ratio (score) identity.
"""
import math
import numpy as np


def _f(g, T, clipped):
    return np.clip(g / (2 * T), -1.0, 1.0) if clipped else np.tanh(g / (2 * T))


def _score(g, s_new, f, T, clipped):
    """d/dg log W(s_new | g) with W = (1 + s_new f(g))/2."""
    if clipped:
        fp = np.where(np.abs(g) < 2 * T, 1.0 / (2 * T), 0.0)
        den = 1.0 + s_new * f
        return np.where(den > 1e-12, s_new * fp / np.maximum(den, 1e-12), 0.0)
    return s_new * (1.0 - s_new * f) / (2 * T)


def solve(T, q, beta=None, lam=None, M=20000, clipped=True, seed=0, dtype=np.float32, estimator='rb'):
    """T, q: raw per-step schedules (length S). beta: raw lag-1 correction per step (TEC: beta=-J_v), or
    lam: Onsager mode, beta_t = lam_t * c_{t-1} with c the self-averaging coefficient P_lin/(2T) (raw: n_lin/(2T)).
    Returns per-time observables of the N -> inf process (N only sets the raw-unit scaling, here N=2000)."""
    N = 2000; sq = math.sqrt(N)
    T = np.asarray(T, float) / sq; q = np.asarray(q, float) / sq
    S = len(T)
    beta = np.zeros(S) if beta is None else np.asarray(beta, float) / sq
    rng = np.random.default_rng(seed)
    Z = np.empty((S + 1, M), dtype)          # iid normals generating eta
    sh = np.empty((S + 1, M), dtype)         # spins s(0..S)
    sc = np.empty((S, M), dtype)             # scores of decisions 0..S-1
    L = np.zeros((S + 1, S + 1))
    Li = np.zeros((S + 1, S + 1))       # inverse of the Cholesky factor, extended one row per step
    G = np.zeros((S + 1, S + 1))
    C = np.zeros((S + 1, S + 1))
    out = dict(su=np.zeros(S + 1), flips=np.zeros(S), plin=np.zeros(S), beta=np.zeros(S), c=np.zeros(S), m_abs=np.zeros(S + 1))
    sh[0] = rng.choice(np.array([-1.0, 1.0], dtype), size=M)
    c_prev = 0.0
    for t in range(S + 1):
        st = sh[t]
        C[t, :t + 1] = sh[:t + 1] @ st / M
        C[t, t] = 1.0
        C[:t, t] = C[t, :t]
        if t == 0:
            L[0, 0] = 1.0; Li[0, 0] = 1.0
        else:
            l = Li[:t, :t] @ C[t, :t]
            L[t, :t] = l
            L[t, t] = math.sqrt(max(1.0 - float(l @ l), 1e-10))
            Li[t, :t] = -(l @ Li[:t, :t]) / L[t, t]
            Li[t, t] = 1.0 / L[t, t]
        Z[t] = rng.standard_normal(M).astype(dtype)
        eta = L[t, :t + 1].astype(dtype) @ Z[:t + 1]
        if t > 0:
            if estimator == 'ibp':   # Gaussian integration by parts: E[s(t) eta(<t)] = C G(t,<t), eta = L Z  =>  G = Li^T E[Z s(t)]
                G[t, :t] = Li[:t, :t].T @ (Z[:t] @ (st - st.mean()) / M)
            elif estimator == 'rb':  # Rao-Blackwellized score identity: E[s(t)|past] = f(g(t-1)); exact lag-1 term E[f'(g(t-1))]
                G[t, t - 1] = c_prev
                if t > 1:
                    G[t, :t - 1] = sc[:t - 1] @ (f_prev - f_prev.mean()) / M
            else:                    # plain likelihood-ratio (score) identity
                G[t, :t] = sc[:t] @ (st - st.mean()) / M
            react = G[t, :t].astype(dtype) @ sh[:t]
        else:
            react = 0.0
        u = eta + react
        out['su'][t] = float(np.mean(st * u))
        if t == S:
            break
        if lam is not None:
            b = (lam[t] * c_prev) if t > 0 else 0.0
        else:
            b = beta[t] if t > 0 else 0.0
        g = u + q[t] * st - (b * sh[t - 1] if t > 0 else 0.0)
        f = _f(g, T[t], clipped)
        s_new = np.where(rng.random(M) < 0.5 * (1.0 + f), 1.0, -1.0).astype(dtype)
        sc[t] = _score(g, s_new, f, T[t], clipped)
        f_prev = f.astype(dtype)
        sh[t + 1] = s_new
        plin = float(np.mean(np.abs(g) < 2 * T[t])) if clipped else float(np.mean(1 - f * f))
        c_prev = plin / (2 * T[t])        # = G(t+1, t) exactly (rescaled); raw c = n_lin/(2T) = sqrt(N) * this
        out['flips'][t] = float(np.mean(s_new != st)); out['plin'][t] = plin
        out['beta'][t] = b * sq; out['c'][t] = c_prev * sq
    out['G'] = G; out['C'] = C
    return out


def cut_from_su(su, N=2000, sumw=-1040.0):
    """cut = (sumw + 1/2 sum_i s_i h_i)/2 with sum_i s_i h_i = N^{3/2} E[s u]."""
    return (sumw + 0.5 * N ** 1.5 * np.asarray(su)) / 2


def simulate(J, T, q, beta=None, lam=None, B=256, clipped=True, seed=0):
    """Finite-N reference with the same parametrization; records the same observables per step."""
    rng = np.random.default_rng(seed); N = len(J); S = len(T)
    beta = np.zeros(S) if beta is None else np.asarray(beta, float)
    s = rng.choice(np.array([-1.0, 1.0], np.float32), size=(B, N)); h = s @ J
    sp = None; c_prev = 0.0
    out = dict(su=np.zeros(S + 1), flips=np.zeros(S), plin=np.zeros(S), c1=np.zeros(S), cuts=None)
    for t in range(S):
        out['su'][t] = float((s * h).sum(1).mean() / N ** 1.5)
        b = 0.0 if sp is None else (lam[t] * c_prev if lam is not None else beta[t])
        g = h + q[t] * s - (b * sp if sp is not None else 0.0)
        f = _f(g, T[t], clipped)
        s_new = np.where(rng.random(s.shape, dtype=np.float32) < 0.5 * (1.0 + f), 1.0, -1.0).astype(np.float32)
        nlin = (np.abs(g) < 2 * T[t]).sum(1) if clipped else (1 - f * f).sum(1)
        c_prev = nlin[:, None] / (2 * T[t])          # per-run online coefficient, as in hardware
        out['flips'][t] = float((s_new != s).mean()); out['plin'][t] = float(nlin.mean() / N)
        out['c1'][t] = float((s_new * s).mean())
        d = s_new - s; sp = s; s = s_new; h = h + d @ J
    out['su'][S] = float((s * h).sum(1).mean() / N ** 1.5)
    out['cuts'] = (-1040.0 + 0.5 * (s * h).sum(1)) / 2
    return out


def _pass(T, q, beta, lam, C_in, G_in, M, clipped, rng, dtype):
    """One sampling pass of the effective process with C, G held fixed (rescaled units). Returns empirical C, G, observables."""
    S = len(T)
    Lc = np.linalg.cholesky(C_in + 1e-9 * np.eye(S + 1)).astype(dtype)
    eta = Lc @ rng.standard_normal((S + 1, M)).astype(dtype)
    Gd = G_in.astype(dtype)
    sh = np.empty((S + 1, M), dtype); F = np.zeros((S + 1, M), dtype); sc = np.zeros((S + 1, M), dtype)
    sh[0] = rng.choice(np.array([-1.0, 1.0], dtype), size=M)
    su = np.zeros(S + 1); flips = np.zeros(S); plin = np.zeros(S); bout = np.zeros(S); cvec = np.zeros(S + 1)
    c_prev = 0.0
    for t in range(S + 1):
        u = eta[t] + (Gd[t, :t] @ sh[:t] if t > 0 else 0.0)
        su[t] = float(np.mean(sh[t] * u))
        if t == S:
            break
        b = (lam[t] * c_prev if lam is not None else beta[t]) if t > 0 else 0.0
        g = u + q[t] * sh[t] - (b * sh[t - 1] if t > 0 else 0.0)
        f = _f(g, T[t], clipped)
        s_new = np.where(rng.random(M) < 0.5 * (1.0 + f), 1.0, -1.0).astype(dtype)
        sc[t] = _score(g, s_new, f, T[t], clipped)
        F[t + 1] = f; sh[t + 1] = s_new
        pl = float(np.mean(np.abs(g) < 2 * T[t])) if clipped else float(np.mean(1 - f * f))
        c_prev = pl / (2 * T[t]); cvec[t + 1] = c_prev
        flips[t] = float(np.mean(s_new != sh[t])); plin[t] = pl; bout[t] = b
    C = (sh @ sh.T).astype(float) / M
    Fc = F - F.mean(1, keepdims=True)
    G = (Fc @ sc.T).astype(float) / M                  # G[t, t'] = E[(f(g(t-1)) - mean) score(t')], valid for t' <= t-2
    G = np.tril(G, -2)
    G[np.arange(1, S + 1), np.arange(0, S)] = cvec[1:]  # exact lag-1 response E[f'(g(t-1))]
    return C, G, dict(su=su, flips=flips, plin=plin, beta=bout, c=cvec[1:])


def solve_iter(T, q, beta=None, lam=None, M=50000, iters=30, alpha=0.3, avg_last=10, clipped=True, seed=0,
               dtype=np.float32, N=2000, init=None):
    """Damped fixed-point iteration of the DMFT map with Polyak averaging of the last iterations."""
    sq = math.sqrt(N)
    T = np.asarray(T, float) / sq; q = np.asarray(q, float) / sq; S = len(T)
    beta = np.zeros(S) if beta is None else np.asarray(beta, float) / sq
    rng = np.random.default_rng(seed)
    if init is None:
        C = np.eye(S + 1); G = np.zeros((S + 1, S + 1))
    else:
        C, G = init
    hist = []; acc = None
    for k in range(iters):
        Cn, Gn, obs = _pass(T, q, beta, lam, C, G, M, clipped, rng, dtype)
        dC = float(np.abs(Cn - C).max()); dG = float(np.abs(Gn - G).max())
        C = (1 - alpha) * C + alpha * Cn; G = (1 - alpha) * G + alpha * Gn
        hist.append(dict(it=k, dC=dC, dG=dG, su_final=float(obs['su'][-1])))
        if k >= iters - avg_last:
            acc = {kk: v.copy() for kk, v in obs.items()} if acc is None else {kk: acc[kk] + obs[kk] for kk in acc}
    out = {kk: v / avg_last for kk, v in acc.items()}
    out['beta'] = out['beta'] * sq; out['c'] = out['c'] * sq
    out['C'] = C; out['G'] = G; out['hist'] = hist
    return out


def _psd_corr(C, floor=1e-8):
    """Nearest-ish PSD correlation matrix (eigenvalue floor, unit diagonal); needed because rows are averaged separately."""
    try:
        return np.linalg.cholesky(C + 1e-9 * np.eye(len(C))) @ np.linalg.cholesky(C + 1e-9 * np.eye(len(C))).T
    except np.linalg.LinAlgError:
        w, V = np.linalg.eigh((C + C.T) / 2)
        X = (V * np.maximum(w, floor)) @ V.T
        d = np.sqrt(np.diag(X))
        return X / d[:, None] / d[None, :] + 1e-9 * np.eye(len(C))


def _pass_upto(T, q, beta, lam, C_in, G_in, M, clipped, rng, dtype, t_end):
    """Sampling pass over times 0..t_end (decisions 0..t_end-1) with C, G fixed. Returns C, G blocks and observables."""
    n = t_end + 1
    Lc = np.linalg.cholesky(_psd_corr(C_in[:n, :n])).astype(dtype)
    eta = Lc @ rng.standard_normal((n, M)).astype(dtype)
    Gd = G_in[:n, :n].astype(dtype)
    sh = np.empty((n, M), dtype); F = np.zeros((n, M), dtype); sc = np.zeros((n, M), dtype)
    sh[0] = rng.choice(np.array([-1.0, 1.0], dtype), size=M)
    su = np.zeros(n); flips = np.zeros(t_end); plin = np.zeros(t_end); bout = np.zeros(t_end); cvec = np.zeros(n)
    c_prev = 0.0
    for t in range(n):
        u = eta[t] + (Gd[t, :t] @ sh[:t] if t > 0 else 0.0)
        su[t] = float(np.mean(sh[t] * u))
        if t == t_end:
            break
        b = (lam[t] * c_prev if lam is not None else beta[t]) if t > 0 else 0.0
        g = u + q[t] * sh[t] - (b * sh[t - 1] if t > 0 else 0.0)
        f = _f(g, T[t], clipped)
        s_new = np.where(rng.random(M) < 0.5 * (1.0 + f), 1.0, -1.0).astype(dtype)
        sc[t] = _score(g, s_new, f, T[t], clipped)
        F[t + 1] = f; sh[t + 1] = s_new
        pl = float(np.mean(np.abs(g) < 2 * T[t])) if clipped else float(np.mean(1 - f * f))
        c_prev = pl / (2 * T[t]); cvec[t + 1] = c_prev
        flips[t] = float(np.mean(s_new != sh[t])); plin[t] = pl; bout[t] = b
    C = (sh @ sh.T).astype(float) / M
    G = ((F - F.mean(1, keepdims=True)) @ sc.T).astype(float) / M
    G = np.tril(G, -2)
    G[np.arange(1, n), np.arange(0, n - 1)] = cvec[1:]
    return C, G, dict(su=su, flips=flips, plin=plin, beta=bout, c=cvec[1:])


def solve_causal(T, q, beta=None, lam=None, M=30000, W=40, it_window=16, warm=4, final_passes=8, clipped=True, seed=0,
                 dtype=np.float32, N=2000):
    """Causal block iteration: rows of C, G for times in [k*W, (k+1)*W] are iterated (damped for `warm` passes, then
    1/n running averages) while all earlier rows keep 1/n running averages. Exploits that the solution up to t does not
    depend on later times. Observables are averaged over `final_passes` full-length passes with the converged C, G."""
    sq = math.sqrt(N)
    T = np.asarray(T, float) / sq; q = np.asarray(q, float) / sq; S = len(T)
    beta = np.zeros(S) if beta is None else np.asarray(beta, float) / sq
    rng = np.random.default_rng(seed)
    C = np.eye(S + 1); G = np.zeros((S + 1, S + 1)); cnt = np.zeros(S + 1)   # per-row averaging counts
    t_lo = 0
    while t_lo < S:
        t_end = min(t_lo + W, S)
        for it in range(it_window):
            Cn, Gn, _ = _pass_upto(T, q, beta, lam, C, G, M, clipped, rng, dtype, t_end)
            for t in range(t_end + 1):
                if t >= t_lo and it < warm:
                    w = 0.5                       # fast damped start for the new window
                else:
                    cnt[t] += 1; w = 1.0 / cnt[t]
                C[t, :t + 1] = (1 - w) * C[t, :t + 1] + w * Cn[t, :t + 1]; C[:t + 1, t] = C[t, :t + 1]
                G[t, :t] = (1 - w) * G[t, :t] + w * Gn[t, :t]
            np.fill_diagonal(C, 1.0)
        t_lo = t_end
    acc = None
    for k in range(final_passes):
        _, _, obs = _pass_upto(T, q, beta, lam, C, G, M, clipped, rng, dtype, S)
        acc = {kk: v.copy() for kk, v in obs.items()} if acc is None else {kk: acc[kk] + obs[kk] for kk in acc}
    out = {kk: v / final_passes for kk, v in acc.items()}
    out['beta'] = out['beta'] * sq; out['c'] = out['c'] * sq; out['C'] = C; out['G'] = G
    return out
