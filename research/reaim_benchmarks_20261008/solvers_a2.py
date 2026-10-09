"""Addendum 2 solver dispatch: solvers.py (frozen) for every method, except ReAIM ASA, which here takes ReAIM's own
options: F in {'max', 'min'} (the |N|-FIFO reduction of Algorithm 2 line 22, ReAIM Step 4) and the initial relative
temperature T0 (ReAIM Table II: GPP 1 -> 0.01, TSP 0.5 -> 0.1, Max-Cut 1 -> 0.1). With F = 'max' and T0 = 1.0,
reaim_a2 performs exactly the operations of solvers.reaim (bit-identical; checked by verify_a2.py)."""
import numpy as np

import solvers as SV


def reaim_a2(P, B, cfg, rng):
    F = {'max': np.max, 'min': np.min}[cfg.get('F', 'max')]
    S = cfg['S']; kset = cfg['kset']; t0, t1 = float(cfg.get('T0', 1.0)), cfg['T1']; it_trial, it_run = 32, 96
    N = P.N; R = len(kset)
    alpha = (t1 / t0) ** (1.0 / (S - 1))
    x = rng.choice(np.array([-1.0, 1.0], dtype=np.float32), size=(B, N))
    cur = SV.Batch(P, x, np.full((B, 20), N, dtype=np.int64), F)
    t, T = 0, t0
    while t < S:
        n = min(it_trial, S - t)
        rep = SV.Batch(P, np.repeat(cur.x, R, 0), np.repeat(cur.fifo, R, 0), F, h=np.repeat(cur.h, R, 0),
                       M=np.repeat(cur.M, R, 0))
        kvec = np.tile(np.array(kset), B)
        for u in range(n):
            rep.step(kvec, T, rng); T *= alpha
        Elast = rep.energy().reshape(B, R)
        pick = Elast.argmin(1)
        t += n
        sel = np.arange(B) * R + pick
        cur = SV.Batch(P, rep.x[sel].copy(), rep.fifo[sel].copy(), F, h=rep.h[sel].copy(), M=rep.M[sel].copy())
        kbest = np.array(kset)[pick]
        n = min(it_run, S - t)
        for u in range(n):
            cur.step(kbest, T, rng); T *= alpha
        t += n
    return cur.x


def run_method(P, cfg, B, seed):
    """solvers.run_method, with ReAIM ASA routed to reaim_a2 (same RNG creation and use)."""
    if cfg['family'] == 'ReAIM':
        rng = np.random.default_rng(seed)
        return reaim_a2(P, B, cfg, rng), None, None
    return SV.run_method(P, cfg, B, seed)
