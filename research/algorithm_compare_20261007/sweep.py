"""PROTOCOL.md sweep: per (method, S) pilot over the method's grid (64 runs), select by pilot mean final cut, then a held-out
final run (256 runs). Usage: python3 sweep.py [S ...]   (default: 1000 250 500 2000 4000). Results: results/S<S>.json."""
import itertools
import json
import math
import sys
import time
from pathlib import Path

import numpy as np
import methods as M

ROOT = Path(__file__).resolve().parent
OUT = ROOT / 'results'
S_LIST = (250, 500, 1000, 2000, 4000)
B_PILOT, B_FINAL = 64, 256


def grids(S):
    P = itertools.product
    return [
        ('SA', [dict(family='SA', T0=a, T1=b, S=S) for a, b in P((20.0, 30.0, 50.0, 80.0), (0.5, 1.0, 2.0, 3.0))]),
        ('SCA', [dict(family='plain', q=q, T0=t, S=S) for q, t in P((2.0, 4.0, 6.0, 8.0), (15.0, 20.0, 30.0, 40.0))]),
        ('TEC', [dict(family='tec', jv=j, q=q, T0=t, S=S) for j, q, t in P((-4.0, -2.0, 2.0, 4.0), (6.0, 8.0), (20.0, 30.0))]),
        ('APC-SCA', [dict(family='apc', q_reset=a, r_q=b, q_lim=c, T0=t, S=S)
                     for a, b, c, t in P((8.0, 16.0, 32.0), (0.9, 0.97), (2.0, 4.0), (15.0, 20.0, 30.0))]),
        ('ReAIM ASA', [dict(family='ReAIM', kset=k, T1=t1, S=S)
                       for k, t1 in P(((64, 128, 256, 512), (128, 256, 512, 1024), (256, 512, 1024, 2000)), (0.1, 0.05))]),
        ('aSB', [dict(family='aSB', dt=d, xi=x, S=S) for d, x in P((0.5, 0.9, 1.25), (0.5, 1.0, 1.5))]),
        ('bSB', [dict(family='bSB', dt=d, xi=x, S=S) for d, x in P((0.5, 0.75, 1.0, 1.25), (0.5, 1.0, 2.0))]),
        ('dSB', [dict(family='dSB', dt=d, xi=x, S=S) for d, x in P((0.5, 0.75, 1.0, 1.25), (0.5, 1.0, 2.0))]),
        ('TEC-T (ours)', [dict(family='tecT', kappa=k, q=q, T0=t, ramp=True, S=S)
                          for k, q, t in P((1.25, 1.75, 2.0), (6.0, 8.0), (12.0, 15.0, 20.0))]),
        ('Onsager SCA (ours)', [dict(family='onsager', lam=l, q=q, T0=t, ramp=True, S=S)
                                for l, q, t in P((0.7, 0.9, 1.05), (6.0, 8.0), (12.0, 15.0, 20.0))]),
    ]


def wilson(k, n, z=1.959963984540054):
    p = k / n; d = 1 + z * z / n
    c = (p + z * z / (2 * n)) / d; h = z * math.sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / d
    return max(0.0, c - h), min(1.0, c + h)


def mcs99(S, p):
    return math.inf if p <= 0 else (S if p >= 1 else S * math.log(0.01) / math.log1p(-p))


def main():
    S_run = [int(a) for a in sys.argv[1:]] or [1000, 250, 500, 2000, 4000]
    np.seterr(all='ignore')   # aSB divergence is detected and flagged invalid, not warned about
    J, sumw = M.m.load()
    OUT.mkdir(exist_ok=True)
    for S in S_run:
        si = S_LIST.index(S); path = OUT / f'S{S}.json'
        res = json.loads(path.read_text()) if path.exists() else {}
        for mi, (name, grid) in enumerate(grids(S)):
            if name in res and 'final' in res[name]:
                continue
            t0 = time.time(); pilot = []
            for gi, cfg in enumerate(grid):
                E, cuts = M.run_method(J, sumw, cfg, B_PILOT, 100000 * mi + 80000 + 1000 * si + gi)
                pilot.append(dict(cfg=cfg, mean_cut=float(cuts.mean()), p33000=float((cuts >= 33000).mean()),
                                  valid=bool(np.isfinite(cuts).all())))
            ok = [i for i, r in enumerate(pilot) if r['valid']]
            sel = max(ok, key=lambda i: (pilot[i]['mean_cut'], -i))
            cfg = grid[sel]
            E, cuts = M.run_method(J, sumw, cfg, B_FINAL, 100000 * mi + 90000 + 1000 * si)
            k = int((cuts >= 33000).sum()); p = k / B_FINAL
            fin = dict(cfg=cfg, runs=B_FINAL, valid=bool(np.isfinite(cuts).all()), mean_cut=float(cuts.mean()),
                       sd_cut=float(cuts.std(ddof=1)), max_cut=float(cuts.max()), k33000=k, p33000=p,
                       p33000_wilson95=wilson(k, B_FINAL), p33200=float((cuts >= 33200).mean()),
                       p33337=float((cuts >= 33337).mean()), mcs99=mcs99(S, p),
                       H_mean=E.mean(1).tolist(), H_p10=np.percentile(E, 10, axis=1).tolist(),
                       H_p90=np.percentile(E, 90, axis=1).tolist(), cuts=cuts.tolist())
            res[name] = dict(method_index=mi, pilot=pilot, selected=sel, final=fin, sec=time.time() - t0)
            path.write_text(json.dumps(res) + '\n')
            print(f"S={S} {name:20s} sel={sel:2d} pilot_mean={pilot[sel]['mean_cut']:.0f} | final mean={fin['mean_cut']:.1f} "
                  f"p33000={p:.3f} p33200={fin['p33200']:.3f} mcs99={fin['mcs99']:.0f} ({time.time() - t0:.0f}s)", flush=True)


if __name__ == '__main__':
    main()
