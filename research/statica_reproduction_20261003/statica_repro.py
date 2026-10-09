"""Reproduce STATICA's K2000 SCA results (JSSC 2021, Table II / Figs. 22, 24, 25). See PROTOCOL.md."""
import hashlib
import json
import math
import sys
import time
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parent
RUD = ROOT.parents[1] / "fpga/v80_snowball/data/WK2000_1.rud"
RUD_SHA256 = "9ed615e5e18726914f12740b7f9bedb6b69676477ed5958fac42c8ab0c252ba7"
ROOT_SEED = 20261003
Q, T_FIN, TARGET_CUT = 4.0, 5.0, 33000
GRID_S = [360, 560, 760, 960, 1160, 1360, 1560]
GRID_TINIT = [50.0, 40.0, 30.0]
TABLE2 = {"short": (560, 30.0), "long": (1560, 40.0)}
PAPER = {"short": dict(mean_H=-66540, mean_cut=32750, p=0.07, t_ms=0.13, tts_ms=8.23),
         "long": dict(mean_H=-67186, mean_cut=33073, p=0.77, t_ms=0.48, tts_ms=1.50)}


def load_graph():
    raw = RUD.read_bytes()
    if hashlib.sha256(raw).hexdigest() != RUD_SHA256:
        raise RuntimeError("WK2000_1.rud hash mismatch")
    lines = raw.split(b"\n", 1)
    n, m = map(int, lines[0].split())
    e = np.array(lines[1].split(), dtype=np.int64).reshape(-1, 3)
    if len(e) != m:
        raise RuntimeError("edge count mismatch")
    i, j, w = e[:, 0] - 1, e[:, 1] - 1, e[:, 2]
    W = np.zeros((n, n), dtype=np.float32)
    W[i, j] = w
    W[j, i] = w
    if np.any(i == j):
        raise RuntimeError("self loop in edge list")
    return -W, (i, j, w)  # J = -w, zero diagonal


def stay_probability(z, T):
    """Eq. (7): clipped-linear non-transition probability."""
    return np.clip(z / (4.0 * T) + 0.5, 0.0, 1.0)


def sca(J, S, t_init, t_fin, q, sigma0, uniforms):
    """Vectorized Algorithm 1 over a batch of independent trials. Returns final spins and flips per step."""
    sigma = sigma0.astype(np.float32).copy()
    h = sigma @ J
    r_T = (t_fin / t_init) ** (1.0 / (S - 1))
    T = t_init
    flips = np.empty((S, len(sigma)), dtype=np.int32)
    for s in range(S):
        flip = stay_probability(h * sigma + q, T) < uniforms(s)
        d = np.where(flip, -2.0 * sigma, 0.0).astype(np.float32)
        sigma += d
        h += d @ J
        flips[s] = flip.sum(axis=1)
        T *= r_T
    return sigma, flips


def sca_scalar(J, S, t_init, t_fin, q, sigma0, U):
    """Literal Algorithm 1 for one trial (lines 1-21), for the equivalence test."""
    N = len(sigma0)
    sigma = [float(x) for x in sigma0]
    r_T = (t_fin / t_init) ** (1.0 / (S - 1))
    T = t_init
    for s in range(S):
        tau = list(sigma)
        for x in range(N):
            hx = sum(float(J[x, y]) * sigma[y] for y in range(N))
            P = min(1.0, max(0.0, (hx * sigma[x] + q) / (4.0 * T) + 0.5))
            if P < U[s, x]:
                tau[x] = -sigma[x]
        sigma = tau
        T *= r_T
    return np.array(sigma)


def self_test():
    rng = np.random.default_rng(7)
    for trial in range(6):
        N, S = 12, 40
        W = np.triu(rng.choice([-1.0, 1.0], size=(N, N)), 1)
        J = (-(W + W.T)).astype(np.float32)
        sigma0 = rng.choice([-1.0, 1.0], size=(1, N))
        U = rng.random((S, 1, N))
        t_init = rng.uniform(2.0, 10.0)
        q = 1.7 + 0.3 * trial
        vec, _ = sca(J, S, t_init, 0.5, q, sigma0, lambda s: U[s])
        ref = sca_scalar(J, S, t_init, 0.5, q, sigma0[0], U[:, 0, :])
        if not np.array_equal(vec[0], ref):
            raise RuntimeError("vectorized SCA differs from literal Algorithm 1")
    return "PASS (6 random 12-spin graphs, 40 steps, shared uniforms)"


def energies(J, sigma):
    return -0.5 * np.einsum("bi,bi->b", sigma, sigma @ J).astype(np.float64)


def cut_from_edges(edges, sigma, chunk=32):
    i, j, w = edges
    out = np.empty(len(sigma), dtype=np.int64)
    for a in range(0, len(sigma), chunk):
        s = sigma[a:a + chunk]
        out[a:a + chunk] = ((s[:, i] != s[:, j]) * w).sum(axis=1)
    return out


def cycles_flips_only(flips, N):
    """Step 1: N/2 cycles for the full field; step s>=2: ceil(F_{s-1}/2) cycles (two flipped spins per cycle)."""
    return N / 2 + np.ceil(flips[:-1] / 2.0).sum(axis=0)


def binom_cdf(k, n, p):
    """P(X <= k) for X ~ Binomial(n, p), summed in log space."""
    if p <= 0:
        return 1.0
    if p >= 1:
        return 0.0 if k < n else 1.0
    lp, lq = math.log(p), math.log1p(-p)
    terms = [math.lgamma(n + 1) - math.lgamma(i + 1) - math.lgamma(n - i + 1) + i * lp + (n - i) * lq
             for i in range(k + 1)]
    m = max(terms)
    return min(1.0, math.exp(m) * sum(math.exp(t - m) for t in terms))


def clopper_pearson(k, n, a=0.05):
    """Exact two-sided interval by bisection on the binomial tails."""
    def solve(f):  # f decreasing in p; find f(p) = 0
        lo, hi = 0.0, 1.0
        for _ in range(200):
            mid = (lo + hi) / 2
            lo, hi = (mid, hi) if f(mid) > 0 else (lo, mid)
        return (lo + hi) / 2
    low = 0.0 if k == 0 else solve(lambda p: a / 2 - (1 - binom_cdf(k - 1, n, p)))
    high = 1.0 if k == n else solve(lambda p: binom_cdf(k, n, p) - a / 2)
    return float(low), float(high)


def tts(t, p):
    return float("nan") if p <= 0 or p >= 1 else t * math.log(0.01) / math.log(1 - p)


def run_config(J, sumw, label, S, t_init, trials, seed_seq, batch=1024):
    sigmas, flips_all = [], []
    batch_seeds = seed_seq.spawn((trials + batch - 1) // batch)
    for b0 in range(0, trials, batch):
        rng = np.random.default_rng(batch_seeds[b0 // batch])
        B = min(batch, trials - b0)
        sigma0 = rng.choice(np.array([-1.0, 1.0], dtype=np.float32), size=(B, J.shape[0]))
        sig, fl = sca(J, S, t_init, T_FIN, Q, sigma0,
                      lambda s: rng.random((B, J.shape[0]), dtype=np.float32))
        sigmas.append(sig.astype(np.int8))
        flips_all.append(fl)
    sigma = np.concatenate(sigmas)
    flips = np.concatenate(flips_all, axis=1)
    H = energies(J, sigma.astype(np.float32))
    cut = (sumw - H) / 2
    k = int((cut >= TARGET_CUT).sum())
    return sigma, flips, dict(label=label, S=S, T_init=t_init, T_fin=T_FIN, q=Q, trials=trials, successes=k,
                              p=k / trials, p_ci95=clopper_pearson(k, trials), mean_H=float(H.mean()),
                              sd_H=float(H.std(ddof=1)), mean_cut=float(cut.mean()), best_cut=float(cut.max()),
                              mean_total_flips=float(flips.sum(axis=0).mean()),
                              mean_cycles_flips_only=float(cycles_flips_only(flips, J.shape[0]).mean()),
                              mean_flips_per_step=flips.mean(axis=1).round(2).tolist())


def main():
    out = {"self_test": self_test()}
    print(out["self_test"], flush=True)
    J, edges = load_graph()
    sumw = float(edges[2].sum())
    out.update(instance=str(RUD.relative_to(ROOT.parents[1])), sha256=RUD_SHA256, sum_w=sumw,
               target_H=sumw - 2 * TARGET_CUT, root_seed=ROOT_SEED)
    root = np.random.SeedSequence(ROOT_SEED)
    children = root.spawn(len(GRID_S) * len(GRID_TINIT) + len(TABLE2))
    results, c = {}, 0
    for name, (S, t_init) in TABLE2.items():
        start = time.perf_counter()
        sigma, flips, r = run_config(J, sumw, name, S, t_init, 4096, children[c]); c += 1
        exact = cut_from_edges(edges, sigma)
        cut = (sumw - energies(J, sigma.astype(np.float32))) / 2
        if not np.array_equal(exact, cut.astype(np.int64)):
            raise RuntimeError(f"edge-list rescoring mismatch in {name}")
        r["edge_rescore"] = f"PASS ({len(sigma)} final states)"
        r["seconds"] = time.perf_counter() - start
        np.save(ROOT / f"final_states_{name}.npy", np.packbits(sigma > 0, axis=1))
        results[name] = r
        print(json.dumps({k: r[k] for k in ("label", "successes", "p", "mean_cut", "seconds")}), flush=True)
    for t_init in GRID_TINIT:
        for S in GRID_S:
            start = time.perf_counter()
            _, _, r = run_config(J, sumw, f"grid_S{S}_T{t_init:g}", S, t_init, 1024, children[c]); c += 1
            r["seconds"] = time.perf_counter() - start
            results[r["label"]] = r
            print(json.dumps({k: r[k] for k in ("label", "successes", "p", "mean_cut", "seconds")}), flush=True)
    out["results"] = results
    (ROOT / "raw_results.json").write_text(json.dumps(out, indent=1) + "\n")


if __name__ == "__main__":
    main()
