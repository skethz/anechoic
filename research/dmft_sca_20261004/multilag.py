"""Track B: multi-lag echo correction, tested with exact finite-N simulation on K2000.

Echo terms: the field of spin i at decision t contains sum_k chat_k(t) s_i(t-k), chat_k(t) = N * dE[s(t)]/dh(t-k) (raw units).
chat_1(t) = n_lin(t-1)/(2T) is the Onsager coefficient used now. chat_k for k >= 2 is measured here by noise injection:
independent fields +-eps on every spin at every step; chat_k(t) = N * E[s_i(t) r_i(t-k)] / eps (linear response).
"""
import json
import sys
import time
from pathlib import Path

import numpy as np

sys.path.insert(0, str(__import__('pathlib').Path(__file__).resolve().parent / '../theory_ideas_20261003'))
import sweeps_abc as m  # noqa: E402

ROOT = Path(__file__).resolve().parent
KMAX = 6


def lam_sched(lam, ramp, S):
    L = np.full(S, float(lam))
    if ramp:
        for t in range(S):
            if t >= 0.7 * S:
                L[t] = lam * (S - 1 - t) / max(1.0, S - 1 - 0.7 * S)
    return L


def run(J, sumw, S, t0, q, lam1, ramp, B, rng, extra=None, eps=0.0):
    """Clipped SCA with online lag-1 Onsager term lam1_t * n_lin(t-1)/(2T) and optional table terms
    extra[k] (array over t) * s(t-k) for k >= 2. eps > 0 injects +-eps fields and returns measured chat_k(t)."""
    N = len(J); T = m.sched(t0, S); L1 = lam_sched(lam1, ramp, S)
    s = rng.choice(np.array([-1.0, 1.0], np.float32), size=(B, N)); h = s @ J
    hist = [s.copy()]                      # s(t), s(t-1), ... up to KMAX
    rbuf = []                              # injected signs r(t-1), r(t-2), ...
    chat = np.zeros((S + 1, KMAX + 1)); c_prev = None; flips = np.zeros(B)
    for t in range(S):
        if eps > 0:
            for k in range(1, min(KMAX, len(rbuf)) + 1):   # control variate: s(t-k) is independent of r(t-k)
                chat[t, k] = N * float(((s - hist[k]) * rbuf[k - 1]).mean()) / eps
        field = h.copy()
        if c_prev is not None and len(hist) > 1:
            field -= L1[t] * c_prev * hist[1]
        if extra is not None:
            for k, tab in extra.items():
                if len(hist) > k:
                    field -= tab[t] * hist[k]
        r = None
        if eps > 0:
            r = rng.choice(np.array([-1.0, 1.0], np.float32), size=(B, N))
            field = field + eps * r
        z = s * field + q
        flip = np.clip(z / (4 * T[t]) + 0.5, 0, 1) < rng.random((B, N), dtype=np.float32)
        c_prev = ((z > -2 * T[t]) & (z < 2 * T[t])).sum(1, keepdims=True) / (2 * T[t])
        d = np.where(flip, -2 * s, 0).astype(np.float32)
        s = s + d; h = h + d @ J; flips += flip.sum(1)
        hist = [s.copy()] + hist[:KMAX]
        if eps > 0:
            rbuf = [r] + rbuf[:KMAX - 1]
    if eps > 0:
        for k in range(1, min(KMAX, len(rbuf)) + 1):
            chat[S, k] = N * float(((s - hist[k]) * rbuf[k - 1]).mean()) / eps
    cut = m.cut(J, sumw, s)
    return cut, flips, chat


def measure(J, sumw, name, S, t0, q, lam1, ramp, B=256, eps=1.0, seed=1):
    rng = np.random.default_rng(seed); t1 = time.time()
    cut, flips, chat = run(J, sumw, S, t0, q, lam1, ramp, B, rng, eps=eps)
    rows = [(t, chat[max(KMAX + 1, t - 20):t, 1:].mean(0)) for t in (S // 8, S // 4, S // 2, 3 * S // 4, 7 * S // 8, S)]
    print(f"{name}: p={float((cut >= 33000).mean()):.3f} ({time.time()-t1:.0f}s); chat_k(t) for k=1..{KMAX}:", flush=True)
    for t, v in rows:
        print(f"   t={t:4d}: " + " ".join(f"{x:8.2f}" for x in v), flush=True)
    return chat


if __name__ == '__main__':
    J, sumw = m.load()
    out = {}
    out['P1'] = measure(J, sumw, 'P1 plain q4 T30 S1560', 1560, 30.0, 4.0, 0.0, False).tolist()
    out['O1'] = measure(J, sumw, 'O1 onsager q8 lam1.05 ramp T12 S960', 960, 12.0, 8.0, 1.05, True).tolist()
    (ROOT / 'multilag_chat.json').write_text(json.dumps(out))
