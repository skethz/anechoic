"""EXPLORATORY, pre-protocol: (1) speed of every method at G-set sizes, (2) where plain SCA's shared parameters
(q, T0, T_fin, all in units of the RMS field sigma) work on SYNTHETIC graphs of each G-set class, (3) the size of the
online Onsager coefficient on sparse graphs. No G-set instance is used. Output: calib.json, calib.log."""
import itertools
import json
import multiprocessing as mp
import sys
import time

import numpy as np

import engine as E
import others as O
import synth

SPECS = {
    'R800+': ('random', 800, 19176, False), 'R800+-': ('random', 800, 19176, True), 'T800+-': ('torus', 20, 40, None),
    'P800+': ('planar', 800, None, False), 'P800+-': ('planar', 800, None, True), 'R2000+': ('random', 2000, 19990, False),
    'T2000+-': ('torus', 40, 50, None), 'P2000+': ('planar', 2000, None, False), 'R1000+': ('random', 1000, 9990, False),
}
_G = {}


def graph(key):
    if key not in _G:
        k, a, b, sg = SPECS[key]
        _G[key] = (synth.random_graph(a, b, sg, 901) if k == 'random' else synth.torus(a, b, 902) if k == 'torus'
                   else synth.planar_like(a, sg, 903))
    return _G[key]


def task(args):
    key, rel, B, seed = args
    g = graph(key); a = g.sigma
    cfg = dict(family='plain', q=rel['q'] * a, T0=rel['T0'] * a, S=rel['S'])
    rng = np.random.default_rng(seed)
    s0 = rng.choice(np.array([-1.0, 1.0], dtype=np.float32), size=(B, g.N))
    t = time.perf_counter(); s, fl = E.run(g, s0, cfg, rng, rel['tfin'] * a); dt = time.perf_counter() - t
    c = g.cut(s)
    return dict(graph=key, **rel, mean_cut=float(c.mean()), max_cut=int(c.max()), mean_flips=float(fl.mean()), sec=dt)


def speed():
    out = {}
    g = graph('R2000+'); B = 64
    rng = np.random.default_rng(1); s0 = rng.choice(np.array([-1.0, 1.0], dtype=np.float32), size=(B, g.N))
    for fam, cfg in (('plain', dict(family='plain', q=8 * g.alpha, T0=20 * g.alpha, S=200)),
                     ('onsager', dict(family='onsager', q=8 * g.alpha, lam=1.0 * g.lam_factor, ramp=True, T0=12 * g.alpha, S=200))):
        E.run(g, s0[:2], dict(cfg, S=5), np.random.default_rng(0), 5 * g.alpha)
        t = time.perf_counter(); E.run(g, s0, cfg, np.random.default_rng(2), 5 * g.alpha); dt = time.perf_counter() - t
        out[fam] = dt / (B * cfg['S'] * g.N) * 1e9
    for fam, fn, cfg in (('SA', O.sa, dict(family='SA', T0=50 * g.alpha, T1=2 * g.alpha, S=200)),
                         ('bSB', O.sb, dict(family='bSB', dt=1.0, xi=1.0, S=200)),
                         ('aSB', O.asb, dict(family='aSB', dt=0.9, xi=1.0, S=200)),
                         ('ReAIM', O.reaim, dict(family='ReAIM', kset=(64, 128, 256, 512), T1=0.1, S=256))):
        fn(g, 2, dict(cfg, S=cfg['S'] if fam == 'ReAIM' else 5), np.random.default_rng(0))
        t = time.perf_counter(); fn(g, B, cfg, np.random.default_rng(3)); dt = time.perf_counter() - t
        out[fam] = dt / (B * cfg['S'] * g.N) * 1e9
    return out   # ns per spin-step (single thread)


def main():
    res = dict(speed_ns_per_spin_step=speed())
    print('speed', json.dumps(res['speed_ns_per_spin_step']), flush=True)
    grid = [dict(q=q, T0=t0, tfin=tf, S=1000) for q, t0, tf in
            itertools.product((0.045, 0.09, 0.18, 0.36, 0.72, 1.44), (0.27, 0.45, 0.9), (0.056, 0.112, 0.224))]
    jobs = [(k, rel, 32, 1000 + i) for k in SPECS for i, rel in enumerate(grid)]
    with mp.Pool(6) as pool:
        rows = []
        for r in pool.imap_unordered(task, jobs, chunksize=2):
            rows.append(r)
    res['plain_scan'] = rows
    for k in SPECS:
        rr = sorted((r for r in rows if r['graph'] == k), key=lambda r: -r['mean_cut'])
        best = rr[0]['mean_cut']; mx = max(r['max_cut'] for r in rr)
        print(f"{k:8s} N={graph(k).N} m={graph(k).m} sigma={graph(k).sigma:.2f} best-found {mx}", flush=True)
        for r in rr[:6]:
            print(f"   q/s {r['q']:.3f} T0/s {r['T0']:.2f} Tf/s {r['tfin']:.3f}: mean {r['mean_cut']:.1f} ({100 * r['mean_cut'] / mx:.2f}% of best-found)"
                  f" flips {r['mean_flips']:.0f}", flush=True)
        # marginal view: best mean cut per q
        for q in sorted({r['q'] for r in rr}):
            m = max(r['mean_cut'] for r in rr if r['q'] == q)
            print(f"      q/s={q:.3f}: best mean {m:.1f} ({100 * m / mx:.2f}%)", flush=True)
    (E.ROOT / 'calib.json').write_text(json.dumps(res, indent=1) + '\n')


if __name__ == '__main__':
    main()
