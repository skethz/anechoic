"""Pilot tuning at the common budget S = 960 steps for the methods without an existing selection (SA, bSB, dSB).
Selection: highest mean final cut over the pilot runs (seed 70001). SCA-family configurations come from the frozen
ablation selections (research/ablation_20261005) and the hardware protocol (O1, P4, T2), i.e. they are not re-tuned here."""
import itertools
import json
import time
from pathlib import Path

import methods as M

ROOT = Path(__file__).resolve().parent
S = 960


def main():
    J, sumw = M.m.load()
    grid = [dict(family='SA', T0=t0, T1=t1, S=S) for t0, t1 in itertools.product((5.0, 10.0, 20.0), (0.25, 0.5, 1.0))]
    grid += [dict(family=f, dt=dt, xi=xi, S=S) for f in ('bSB', 'dSB') for dt, xi in itertools.product((0.5, 0.75, 1.0, 1.25), (0.5, 1.0, 2.0))]
    out = []
    for i, cfg in enumerate(grid):
        B = 32 if cfg['family'] == 'SA' else 128
        t = time.time(); E, cuts = M.run_method(J, sumw, cfg, B, 70001 + i)
        r = dict(cfg, B=B, mean_cut=float(cuts.mean()), max_cut=float(cuts.max()), p33000=float((cuts >= 33000).mean()), sec=time.time() - t)
        out.append(r); print(json.dumps(r), flush=True)
    (ROOT / 'tune.json').write_text(json.dumps(out, indent=1) + '\n')
    for fam in ('SA', 'bSB', 'dSB'):
        best = max((r for r in out if r['family'] == fam), key=lambda r: r['mean_cut'])
        print('BEST', json.dumps(best))


if __name__ == '__main__':
    main()
