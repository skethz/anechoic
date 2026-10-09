"""EXPLORATORY quick tests of two theory-motivated ideas on WK2000_1 (laptop CPU; not a frozen protocol).

A. Onsager-corrected SCA (dynamical TAP for synchronous updates).
   Spin i's value at t-1 echoes back into its field at t through every j:
   echo_i ~ s_i(t-1) * sum_j J_ij^2 chi_j, with chi_j = 1/(2T) inside Eq. 7's linear band and 0 outside.
   For +-1 couplings this is s_i(t-1) * n_lin(t-1) / (2 T_{t-1}).
   The decision uses z = s_i * (h_i - lam * echo_i) + q. lam=0 is plain SCA; lam=1 is the full Onsager subtraction.
B. IAMP stopped at t* <= 1, then a short low-temperature SCA started from sign(m).
"""
import json
import math
import sys
import time
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parent
PROJ = ROOT.parents[1]
sys.path.insert(0, str(PROJ / "research/statica_reproduction_20261003"))
sys.path.insert(0, str(PROJ / "research/iamp_evaluation_20261003"))
import statica_repro as sr
import iamp_eval as ie

Jf, edges = sr.load_graph()
J = Jf.astype(np.float32); N = len(J); SUMW = float(edges[2].sum())
Q = 4.0


def sca(sigma0, S, t0, t1, rng, lam=0.0):
    """Eq. 7 SCA (uniform thresholds) with optional Onsager echo subtraction. Returns final spins, total flips."""
    s = sigma0.astype(np.float32).copy(); h = s @ J
    r = (t1 / t0) ** (1.0 / max(1, S - 1)); T = t0
    s_prev = None; c_prev = None; flips = np.zeros(len(s))
    for _ in range(S):
        field = h if (lam == 0.0 or s_prev is None) else h - lam * c_prev[:, None] * s_prev
        z = s * field + Q
        p_stay = np.clip(z / (4 * T) + 0.5, 0, 1)
        n_lin = ((z > -2 * T) & (z < 2 * T)).sum(1)
        flip = p_stay < rng.random(s.shape, dtype=np.float32)
        d = np.where(flip, -2 * s, 0).astype(np.float32)
        s_prev, c_prev = s.copy(), n_lin / (2 * T)
        s += d; h += d @ J; flips += flip.sum(1); T *= r
    return s, flips


def cut(s):
    return (SUMW + 0.5 * np.einsum("bi,bi->b", s, s @ J)) / 2


def iamp_stop(gauss, K, kstop, table):
    """El Alaoui-Montanari IAMP (same recurrence as iamp_eval), stopped after kstop of K steps."""
    X, D1, D2, G = table; A = J.astype(np.float64) / math.sqrt(N)
    B = len(gauss); dt = 1 / K
    dz = gauss * np.sqrt(N * dt / np.sum(gauss * gauss, 1, keepdims=True))
    m = dz.copy(); x = dz.copy(); z = dz.copy(); prev = np.zeros_like(m); corr = np.zeros_like(m); ap = np.ones(B); hist = [dz.copy()]
    for k in range(1, kstop):
        row = min(255, int(k * 256 / K))
        u = np.interp(x, X, D2[row]); v = G[row] * np.interp(x, X, D1[row])
        u /= np.sqrt(np.mean(u * u, 1, keepdims=True)); am = np.mean(u, 1)
        dz = m @ A + corr - ap[:, None] * prev - z
        for hh in hist[::-1][:5]:
            dz -= hh * np.sum(dz * hh, 1, keepdims=True) / np.sum(hh * hh, 1, keepdims=True)
        dz *= np.sqrt(N * dt / np.sum(dz * dz, 1, keepdims=True))
        mn = m + u * dz; mn *= np.sqrt(N * (k + 1) * dt / np.sum(mn * mn, 1, keepdims=True))
        corr += (am - ap)[:, None] * prev; prev = m; m = mn; ap = am
        x += v * dt + dz; z += dz; hist.append(dz.copy()); hist = hist[-5:]
    return m


def summarize(c, flips, extra=None):
    k = int((c >= 33000).sum())
    d = dict(trials=len(c), successes=k, p=k / len(c), p_ci95=sr.clopper_pearson(k, len(c)),
             mean_cut=float(c.mean()), mean_flips=float(np.mean(flips)))
    d.update(extra or {})
    return d


def main():
    out = {"A": {}, "B": {}}
    # A: Onsager-corrected SCA at both STATICA settings, paired initial states across lam.
    for name, (S, t0) in {"long": (1560, 40.0), "short": (560, 30.0)}.items():
        init = np.random.default_rng(11).choice(np.array([-1.0, 1.0], dtype=np.float32), size=(1024, N))
        for lam in (0.0, 0.5, 1.0, -0.5):
            t = time.perf_counter()
            s, fl = sca(init, S, t0, 5.0, np.random.default_rng(int(1000 + 100 * lam)), lam)
            out["A"][f"{name}_lam{lam}"] = summarize(cut(s), fl, dict(S=S, T0=t0, lam=lam, sec=time.perf_counter() - t))
            print("A", name, lam, json.dumps({k: out['A'][f'{name}_lam{lam}'][k] for k in ('p', 'mean_cut', 'mean_flips')}), flush=True)
    # B: IAMP(t*) + SCA finish
    table = ie.load_table()
    for K, frac in ((64, 1.0), (128, 0.9), (128, 1.0), (256, 0.9)):
        g = np.random.default_rng(K * 10 + int(frac * 10)).normal(size=(512, N))
        t = time.perf_counter(); m = iamp_stop(g, K, int(round(frac * K)), table); ti = time.perf_counter() - t
        s0 = np.where(m >= 0, 1.0, -1.0).astype(np.float32)
        out["B"][f"K{K}_t{frac}_sign"] = summarize(cut(s0), np.zeros(len(s0)), dict(iamp_mvm=int(round(frac * K))))
        for Ta, Sf in ((8.0, 200), (8.0, 400), (12.0, 400)):
            s, fl = sca(s0, Sf, Ta, 5.0, np.random.default_rng(K + Sf + int(Ta)))
            key = f"K{K}_t{frac}_sca_T{Ta}_S{Sf}"
            out["B"][key] = summarize(cut(s), fl, dict(iamp_mvm=int(round(frac * K)), sca_steps=Sf, Ta=Ta))
            print("B", key, json.dumps({k: out['B'][key][k] for k in ('p', 'mean_cut', 'mean_flips')}), flush=True)
    (ROOT / "quick_tests.json").write_text(json.dumps(out, indent=1) + "\n")


if __name__ == "__main__":
    main()
