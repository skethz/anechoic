"""GPU (PyTorch) port of dmft.solve_causal: same effective process, estimators and causal block iteration, larger M."""
import math
import numpy as np
import torch

torch.backends.cuda.matmul.allow_tf32 = False
torch.backends.cudnn.allow_tf32 = False
DEV = torch.device('cuda')


def psd_corr(C, floor=1e-8):
    try:
        np.linalg.cholesky(C + 1e-9 * np.eye(len(C)))
        return C + 1e-9 * np.eye(len(C))
    except np.linalg.LinAlgError:
        w, V = np.linalg.eigh((C + C.T) / 2)
        X = (V * np.maximum(w, floor)) @ V.T
        d = np.sqrt(np.diag(X))
        return X / d[:, None] / d[None, :] + 1e-9 * np.eye(len(C))


def pass_upto(T, q, beta, lam, C_in, G_in, M, clipped, gen, t_end):
    n = t_end + 1
    Lc = torch.linalg.cholesky(torch.from_numpy(psd_corr(C_in[:n, :n])).to(DEV, torch.float64)).float()
    Gd = torch.from_numpy(G_in[:n, :n]).float().to(DEV)
    Z = torch.randn(n, M, device=DEV, generator=gen)
    sh = torch.empty(n, M, device=DEV); F = torch.zeros(n, M, device=DEV); sc = torch.zeros(n, M, device=DEV)
    sh[0] = torch.where(torch.rand(M, device=DEV, generator=gen) < 0.5, 1.0, -1.0)
    su = np.zeros(n); flips = np.zeros(t_end); plin = np.zeros(t_end); cvec = np.zeros(n)
    c_prev = 0.0
    for t in range(n):
        u = Lc[t, :t + 1] @ Z[:t + 1]
        if t > 0:
            u = u + Gd[t, :t] @ sh[:t]
        su[t] = float((sh[t] * u).mean())
        if t == t_end:
            break
        b = (lam[t] * c_prev if lam is not None else beta[t]) if t > 0 else 0.0
        g = u + q[t] * sh[t]
        if t > 0:
            g = g - b * sh[t - 1]
        Tt = T[t]
        if clipped:
            f = torch.clamp(g / (2 * Tt), -1.0, 1.0)
        else:
            f = torch.tanh(g / (2 * Tt))
        s_new = torch.where(torch.rand(M, device=DEV, generator=gen) < 0.5 * (1.0 + f), 1.0, -1.0)
        if clipped:
            fp = (g.abs() < 2 * Tt).float() / (2 * Tt)
            den = 1.0 + s_new * f
            sc[t] = torch.where(den > 1e-12, s_new * fp / den.clamp_min(1e-12), torch.zeros_like(den))
            pl = float((g.abs() < 2 * Tt).float().mean())
        else:
            sc[t] = s_new * (1.0 - s_new * f) / (2 * Tt)
            pl = float((1 - f * f).mean())
        F[t + 1] = f; sh[t + 1] = s_new
        c_prev = pl / (2 * Tt); cvec[t + 1] = c_prev
        flips[t] = float((s_new != sh[t]).float().mean()); plin[t] = pl
    C = (sh @ sh.T).double().cpu().numpy() / M
    G = ((F - F.mean(1, keepdim=True)) @ sc.T).double().cpu().numpy() / M
    del Z, sh, F, sc
    G = np.tril(G, -2)
    G[np.arange(1, n), np.arange(0, n - 1)] = cvec[1:]
    return C, G, dict(su=su, flips=flips, plin=plin, c=cvec[1:])


def solve_causal(T, q, beta=None, lam=None, M=1_000_000, W=40, it_window=32, warm=4, final_passes=4, clipped=True,
                 seed=0, N=2000, log=None):
    sq = math.sqrt(N)
    T = np.asarray(T, float) / sq; q = np.asarray(q, float) / sq; S = len(T)
    beta = np.zeros(S) if beta is None else np.asarray(beta, float) / sq
    gen = torch.Generator(device=DEV); gen.manual_seed(seed)
    C = np.eye(S + 1); G = np.zeros((S + 1, S + 1)); cnt = np.zeros(S + 1)
    t_lo = 0
    while t_lo < S:
        t_end = min(t_lo + W, S)
        for it in range(it_window):
            Cn, Gn, obs = pass_upto(T, q, beta, lam, C, G, M, clipped, gen, t_end)
            for t in range(t_end + 1):
                if t >= t_lo and it < warm:
                    w = 0.5
                else:
                    cnt[t] += 1; w = 1.0 / cnt[t]
                C[t, :t + 1] = (1 - w) * C[t, :t + 1] + w * Cn[t, :t + 1]; C[:t + 1, t] = C[t, :t + 1]
                G[t, :t] = (1 - w) * G[t, :t] + w * Gn[t, :t]
            np.fill_diagonal(C, 1.0)
        if log:
            log(f"window to t={t_end}: su[t_end]={obs['su'][-1]:.4f}")
        t_lo = t_end
    acc = None
    for k in range(final_passes):
        _, _, obs = pass_upto(T, q, beta, lam, C, G, M, clipped, gen, S)
        acc = {kk: v.copy() for kk, v in obs.items()} if acc is None else {kk: acc[kk] + obs[kk] for kk in acc}
    out = {kk: v / final_passes for kk, v in acc.items()}
    out['c'] = out['c'] * sq; out['C'] = C; out['G'] = G
    return out


def solve_causal2(T, q, beta=None, lam=None, M=1_000_000, W=40, it_window=32, warm=6, final_passes=4, clipped=True,
                  seed=0, N=2000, log=None):
    """Like solve_causal, but after warm-up each window replaces the whole leading block of C and G by the running mean of
    that window's pass estimates (common weights), so C stays an exact average of Gram matrices (consistent and PSD)."""
    sq = math.sqrt(N)
    T = np.asarray(T, float) / sq; q = np.asarray(q, float) / sq; S = len(T)
    beta = np.zeros(S) if beta is None else np.asarray(beta, float) / sq
    gen = torch.Generator(device=DEV); gen.manual_seed(seed)
    C = np.eye(S + 1); G = np.zeros((S + 1, S + 1))
    t_lo = 0
    while t_lo < S:
        t_end = min(t_lo + W, S); n = t_end + 1
        Cacc = np.zeros((n, n)); Gacc = np.zeros((n, n)); k = 0
        for it in range(it_window):
            Cn, Gn, obs = pass_upto(T, q, beta, lam, C, G, M, clipped, gen, t_end)
            if it < warm:
                C[t_lo:n, :n] = 0.5 * C[t_lo:n, :n] + 0.5 * Cn[t_lo:n, :n]; C[:n, t_lo:n] = C[t_lo:n, :n].T
                G[t_lo:n, :n] = 0.5 * G[t_lo:n, :n] + 0.5 * Gn[t_lo:n, :n]
                np.fill_diagonal(C, 1.0)
            else:
                k += 1; Cacc += Cn; Gacc += Gn
                C[:n, :n] = Cacc / k; G[:n, :n] = Gacc / k
        if log:
            log(f"window to t={t_end}: su[t_end]={obs['su'][-1]:.4f}")
        t_lo = t_end
    acc = None
    for k in range(final_passes):
        _, _, obs = pass_upto(T, q, beta, lam, C, G, M, clipped, gen, S)
        acc = {kk: v.copy() for kk, v in obs.items()} if acc is None else {kk: acc[kk] + obs[kk] for kk in acc}
    out = {kk: v / final_passes for kk, v in acc.items()}
    out['c'] = out['c'] * sq; out['C'] = C; out['G'] = G
    return out
