"""EXPLORATORY, pre-protocol, SYNTHETIC graphs only: (a) plain-SCA T_fin for the +1 classes (calib2 optimum at the
T_fin edge), (b) whether the K2000-scaled SA grid and the bSB grid are edge-limited on sparse graphs. calib3.json/.log"""
import itertools
import json
import multiprocessing as mp

import numpy as np

import calib
import calib2
import engine as E
import others as O


def task(args):
    key, cfg, B, seed = args
    g = calib.graph(key); a = g.sigma; rng = np.random.default_rng(seed)
    if cfg['family'] == 'plain':
        c = dict(family='plain', q=cfg['q'] * a, T0=cfg['T0'] * a, S=cfg['S'])
        s0 = rng.choice(np.array([-1.0, 1.0], dtype=np.float32), size=(B, g.N))
        s, _ = E.run(g, s0, c, rng, cfg['tfin'] * a)
    elif cfg['family'] == 'SA':
        s = O.sa(g, B, dict(family='SA', T0=cfg['T0'] * a, T1=cfg['T1'] * a, S=cfg['S']), rng)
    else:
        s = O.sb(g, B, cfg, rng)
    cut = g.cut(s)
    return dict(graph=key, cfg=cfg, mean_cut=float(cut.mean()), max_cut=int(cut.max()))


def main():
    S, B = 1000, 32
    jobs = []
    for k in calib2.PLUS:
        for i, (q, t0, tf) in enumerate(itertools.product((1.08, 1.44, 2.16, 2.88), (0.9, 1.35, 1.8, 2.7), (0.224, 0.448, 0.672, 0.896))):
            jobs.append((k, dict(family='plain', q=q, T0=t0, tfin=tf, S=S), B, 4000 + i))
    sa_T0 = tuple(round(x / 44.721, 4) for x in (20, 30, 50, 80)); sa_T1 = tuple(round(x / 44.721, 4) for x in (0.5, 1, 2, 3))
    for k in calib.SPECS:
        for i, (t0, t1) in enumerate(itertools.product(sa_T0, sa_T1)):
            jobs.append((k, dict(family='SA', T0=t0, T1=t1, S=S), B, 5000 + i))
        for i, (dt, xi) in enumerate(itertools.product((0.5, 0.75, 1.0, 1.25), (0.5, 0.75, 1.0, 2.0))):
            jobs.append((k, dict(family='bSB', dt=dt, xi=xi, S=S), B, 6000 + i))
    with mp.Pool(6) as pool:
        rows = list(pool.imap_unordered(task, jobs, chunksize=2))
    for k in calib2.PLUS:
        rr = sorted((r for r in rows if r['graph'] == k and r['cfg']['family'] == 'plain'), key=lambda r: -r['mean_cut'])
        mx = max(r['max_cut'] for r in rr)
        print(f"{k:8s} plain best-found {mx}: " + " | ".join(
            f"q {r['cfg']['q']} T0 {r['cfg']['T0']} Tf {r['cfg']['tfin']}: {r['mean_cut']:.1f}" for r in rr[:4]), flush=True)
        for tf in (0.224, 0.448, 0.672, 0.896):
            print(f"      Tf/s={tf}: best mean {max(r['mean_cut'] for r in rr if r['cfg']['tfin'] == tf):.1f}", flush=True)
    for fam in ('SA', 'bSB'):
        for k in calib.SPECS:
            rr = sorted((r for r in rows if r['graph'] == k and r['cfg']['family'] == fam), key=lambda r: -r['mean_cut'])
            c = rr[0]['cfg']
            print(f"{fam:4s} {k:8s} best {rr[0]['mean_cut']:.1f} at " + (f"T0/s {c['T0']} T1/s {c['T1']}" if fam == 'SA' else f"dt {c['dt']} xi {c['xi']}")
                  + f" | 2nd {rr[1]['mean_cut']:.1f}", flush=True)
    (E.ROOT / 'calib3.json').write_text(json.dumps(rows, indent=1) + '\n')


if __name__ == '__main__':
    main()
