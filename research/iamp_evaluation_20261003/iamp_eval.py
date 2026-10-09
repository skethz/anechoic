"""IAMP K-sweep on WK2000_1 against the reproduced STATICA reference. See PROTOCOL.md."""
import json
import math
import sys
import time
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parent
PROJ = ROOT.parents[1]
sys.path.insert(0, str(PROJ / "research/statica_reproduction_20261003"))
import statica_repro as sr  # graph loader, edge rescoring, Clopper-Pearson

TABLE = PROJ / "research/theory_alternatives_20260929/trial/tables/parisi_zero.npz"
K_GRID = [64, 128, 256, 512, 1024, 2048]
TRIALS, BATCH, ROOT_SEED, MEMORY = 1024, 256, 20261005, 5


def load_table():
    p = np.load(TABLE)
    return p["x"], p["d1"], p["d2"], p["gamma"]


def reference_recurrence(A, gauss, K, table):
    """Literal transcription of the validated CPU reference in continuous.py (numerical_validation), any K."""
    X, D1, D2, G = table
    N = A.shape[0]; dt = 1 / K
    dz = gauss * np.sqrt(N * dt / np.sum(gauss * gauss, axis=1, keepdims=True))
    m = dz.copy(); x = dz.copy(); z = dz.copy(); prev = np.zeros_like(m); corr = np.zeros_like(m)
    ap = np.ones(len(m)); history = [dz.copy()]
    for k in range(1, K):
        row = min(255, int(k * 256 / K))
        u = np.array([np.interp(a, X, D2[row]) for a in x])
        v = G[row] * np.array([np.interp(a, X, D1[row]) for a in x])
        u /= np.sqrt(np.mean(u * u, axis=1, keepdims=True)); amean = np.mean(u, axis=1)
        dz = m @ A + corr - ap[:, None] * prev - z
        for h in history[::-1][:MEMORY]:
            dz -= h * np.sum(dz * h, axis=1, keepdims=True) / np.sum(h * h, axis=1, keepdims=True)
        dz *= np.sqrt(N * dt / np.sum(dz * dz, axis=1, keepdims=True))
        mn = m + u * dz; mn *= np.sqrt(N * (k + 1) * dt / np.sum(mn * mn, axis=1, keepdims=True))
        corr += (amean - ap)[:, None] * prev; prev = m; m = mn; ap = amean
        x += v * dt + dz; z += dz; history.append(dz.copy())
    return m


def iamp(A, J, gauss, K, table, sumw=None):
    """Vectorized IAMP; returns final continuous m, final-state cut, best-visited cut (if sumw given)."""
    X, D1, D2, G = table
    B, N = gauss.shape; dt = 1 / K
    dz = gauss * np.sqrt(N * dt / np.sum(gauss * gauss, axis=1, keepdims=True))
    m = dz.copy(); x = dz.copy(); z = dz.copy(); prev = np.zeros_like(m); corr = np.zeros_like(m)
    ap = np.ones(B); history = [dz.copy()]
    best = None

    def cut_of(mm):
        s = np.where(mm >= 0, 1.0, -1.0)
        return (sumw + 0.5 * np.einsum("bi,bi->b", s, s @ J)) / 2

    if sumw is not None:
        best = cut_of(m)
    for k in range(1, K):
        row = min(255, int(k * 256 / K))
        u = np.interp(x, X, D2[row]); v = G[row] * np.interp(x, X, D1[row])
        u /= np.sqrt(np.mean(u * u, axis=1, keepdims=True)); amean = np.mean(u, axis=1)
        dz = m @ A + corr - ap[:, None] * prev - z
        for h in history[::-1][:MEMORY]:
            dz -= h * np.sum(dz * h, axis=1, keepdims=True) / np.sum(h * h, axis=1, keepdims=True)
        dz *= np.sqrt(N * dt / np.sum(dz * dz, axis=1, keepdims=True))
        mn = m + u * dz; mn *= np.sqrt(N * (k + 1) * dt / np.sum(mn * mn, axis=1, keepdims=True))
        corr += (amean - ap)[:, None] * prev; prev = m; m = mn; ap = amean
        x += v * dt + dz; z += dz; history.append(dz.copy())
        if len(history) > MEMORY:
            history.pop(0)
        if sumw is not None:
            best = np.maximum(best, cut_of(m))
    final = cut_of(m) if sumw is not None else None
    return m, final, best


def self_test(table):
    rng = np.random.default_rng(933)
    worst = 0.0
    for n, K in [(32, 16), (40, 23), (24, 64)]:
        j = np.triu(rng.choice([-1.0, 1.0], size=(n, n)), 1); j = j + j.T
        A = j / math.sqrt(n); g = rng.normal(size=(3, n))
        a = reference_recurrence(A, g.copy(), K, table)
        b, _, _ = iamp(A, j, g.copy(), K, table)
        worst = max(worst, float(np.max(np.abs(a - b))))
    if worst > 1e-10:
        raise RuntimeError(("vectorized IAMP differs from reference recurrence", worst))
    return f"PASS (max |diff| {worst:.1e} over 3 graphs, K up to 64)"


def main():
    table = load_table()
    out = {"self_test": self_test(table)}
    print(out["self_test"], flush=True)
    Jf, edges = sr.load_graph()
    J = Jf.astype(np.float64); N = len(J); A = J / math.sqrt(N)
    sumw = float(edges[2].sum())
    seeds = np.random.SeedSequence(ROOT_SEED).spawn(len(K_GRID))
    res = {}
    for K, ss in zip(K_GRID, seeds):
        start = time.perf_counter()
        finals, bests, rescored = [], [], 0
        for bseed in ss.spawn(TRIALS // BATCH):
            rng = np.random.default_rng(bseed)
            m, final, best = iamp(A, J, rng.normal(size=(BATCH, N)), K, table, sumw)
            if rescored < 256:  # independent edge-list rescoring of final states
                s = np.where(m >= 0, 1, -1).astype(np.int8)
                if not np.array_equal(sr.cut_from_edges(edges, s), final.astype(np.int64)):
                    raise RuntimeError("edge rescoring mismatch")
                rescored += len(s)
            finals.append(final); bests.append(best)
        final = np.concatenate(finals); best = np.concatenate(bests)
        kf, kb = int((final >= sr.TARGET_CUT).sum()), int((best >= sr.TARGET_CUT).sum())
        res[K] = dict(K=K, trials=TRIALS, final_successes=kf, final_p=kf / TRIALS, final_p_ci95=sr.clopper_pearson(kf, TRIALS),
                      final_mean_cut=float(final.mean()), final_sd_cut=float(final.std(ddof=1)),
                      best_successes=kb, best_p=kb / TRIALS, best_p_ci95=sr.clopper_pearson(kb, TRIALS),
                      best_mean_cut=float(best.mean()), best_max_cut=float(best.max()),
                      final_mean_energy_normalized=float(-(sumw - 2 * final.mean()) / N ** 1.5),
                      edge_rescored_final_states=rescored, seconds=time.perf_counter() - start)
        print(json.dumps(res[K]), flush=True)
    out["results"] = res
    (ROOT / "raw_results.json").write_text(json.dumps(out, indent=1) + "\n")


if __name__ == "__main__":
    main()
