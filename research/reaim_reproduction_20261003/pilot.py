"""PILOT, noise-free ReAIM ASA (ISCA 2024, Algorithms 2-3) on WK2000_1. Exploratory, not a frozen reproduction.

Specified by the paper and implemented as written:
- Alg. 2: Delta_i = energy change of flipping i. N = {Delta<0}.
  If |N| > q: th = min positive Delta (greedy PM). Otherwise th = Delta_min + r*T*(Delta_max - Delta_min) (SA-like MN).
  C = {Delta <= th}. Flip min(k, |C|) spins drawn uniformly from C, all at once.
  T <- alpha*T. q <- max of |N| over the past 20 iterations (F = max for Max-Cut, held in a 20-entry FIFO).
- Alg. 3: trial phase runs candidate k values as parallel replicas for ITER_trial iterations. The k with the lowest
  energy continues for ITER_run iterations. H_best is updated after each run phase.
- Max-Cut settings (Table II): T 1 -> 0.1, Ising form.

Not specified for K2000, so these are DECLARED ASSUMPTIONS:
- Candidate k sets, ITER_trial = 32, ITER_run = 96, total iterations swept.
- FIFO initialised with N, so q starts large. If no Delta > 0, th = Delta_max.
- T advances along the continuing trajectory, and replicas inherit the parent's FIFO.
- No ReRAM read or write noise; noise-free is an upper bound on algorithm quality.
"""
import json
import sys
import time
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT.parent / "statica_reproduction_20261003"))
import statica_repro as sr

TARGET = 33000


class Batch:
    """B independent ASA searches; trial replicas expand the batch to B*len(kset)."""

    def __init__(self, J, x, fifo, F=np.max):
        self.J, self.x, self.h, self.fifo, self.F = J, x, x @ J, fifo, F

    def energy(self):
        return -0.5 * np.einsum("bi,bi->b", self.x, self.h)

    def step(self, k, T, rng):
        x, h, J = self.x, self.h, self.J
        d = 2.0 * x * h                                   # Delta_i: H(x^i) - H(x), with H = -x^T J x / 2
        nN = (d < 0).sum(1)
        q = self.F(self.fifo, axis=1)
        pos = np.where(d > 0, d, np.inf).min(1)
        dmin, dmax = d.min(1), d.max(1)
        th_pm = np.where(np.isfinite(pos), pos, dmax)
        th_mn = dmin + rng.random(len(x)) * T * (dmax - dmin)
        th = np.where(nN > q, th_pm, th_mn)
        C = d <= th[:, None]
        keys = np.where(C, rng.random(C.shape), -1.0)
        kk = np.broadcast_to(np.asarray(k), (len(x),))
        kmax = int(kk.max())
        idx = np.argpartition(-keys, kmax - 1, axis=1)[:, :kmax]
        sel = np.take_along_axis(keys, idx, 1) >= 0
        sel &= np.arange(kmax)[None, :] < kk[:, None]      # at most k_b flips for replica b
        rows = np.nonzero(sel)
        flat = idx[rows]
        delta = np.zeros_like(x)
        delta[rows[0], flat] = -2.0 * x[rows[0], flat]
        x += delta
        h += delta @ J
        self.fifo = np.concatenate([self.fifo[:, 1:], nN[:, None]], axis=1)


def asa(J, B, iters, kset, it_trial, it_run, rng, t0=1.0, t1=0.1, F=np.max):
    N = J.shape[0]
    alpha = (t1 / t0) ** (1.0 / (iters - 1))
    x = rng.choice(np.array([-1.0, 1.0], dtype=np.float32), size=(B, N))
    cur = Batch(J, x, np.full((B, 20), N, dtype=np.int64), F)
    best = cur.energy()
    t, T = 0, t0
    R = len(kset)
    while t < iters:
        # Trial phase: replicate each search once per candidate k.
        n = min(it_trial, iters - t)
        rep = Batch(J, np.repeat(cur.x, R, 0), np.repeat(cur.fifo, R, 0), F)
        kvec = np.tile(np.array(kset), B)
        for _ in range(n):
            rep.step(kvec, T, rng); T *= alpha
        t += n
        e = rep.energy().reshape(B, R)
        pick = e.argmin(1)
        sel = np.arange(B) * R + pick
        cur = Batch(J, rep.x[sel].copy(), rep.fifo[sel].copy(), F)
        kbest = np.array(kset)[pick]
        # Run phase with each search's chosen k.
        n = min(it_run, iters - t)
        for _ in range(n):
            cur.step(kbest, T, rng); T *= alpha
        t += n
        best = np.minimum(best, cur.energy())
    return best, cur.energy()


def main():
    Jf, edges = sr.load_graph()
    J = Jf.astype(np.float32)
    sumw = float(edges[2].sum())
    out = {}
    B = 256
    for kset in ((1, 2, 4, 8), (4, 8, 16, 32)):
        for iters in (256, 512, 1024, 2048, 4096):
            rng = np.random.default_rng(hash((kset, iters)) % 2**32)
            t = time.perf_counter()
            best, final = asa(J, B, iters, kset, 32, 96, rng)
            cut = (sumw - best) / 2
            k = int((cut >= TARGET).sum())
            key = f"k{'-'.join(map(str, kset))}_it{iters}"
            out[key] = dict(kset=kset, iters=iters, trials=B, successes=k, p=k / B, p_ci95=sr.clopper_pearson(k, B),
                            mean_best_cut=float(cut.mean()), seconds=time.perf_counter() - t)
            print(key, json.dumps(out[key]), flush=True)
    (ROOT / "pilot_results.json").write_text(json.dumps(out, indent=1) + "\n")


if __name__ == "__main__":
    main()
