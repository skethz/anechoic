"""Final algorithm-level comparison at the common budget S = 960 steps: 256 independent runs per method (seed 71001+i).
Configurations: SA, bSB, dSB tuned here (tune.py, tune_sa_ext.json; selection by mean final cut); SCA family from the
frozen selections used on hardware (P4 plain, T2 TEC, O1 Onsager) and the ablations (APC E1, TEC-T E1)."""
import json
import time
from pathlib import Path

import numpy as np
import methods as M

ROOT = Path(__file__).resolve().parent
S = 960
METHODS = [
    ('SA', dict(family='SA', T0=50.0, T1=3.0, S=S)),
    ('SCA (plain)', dict(family='plain', q=8.0, T0=30.0, S=S)),
    ('TEC', dict(family='tec', q=8.0, jv=-4.0, T0=30.0, S=S)),
    ('APC-SCA', dict(family='apc', q_reset=16.0, r_q=0.9, q_lim=4.0, T0=30.0, S=S)),
    ('dSB', dict(family='dSB', dt=1.0, xi=1.0, S=S)),
    ('bSB', dict(family='bSB', dt=1.0, xi=1.0, S=S)),
    ('TEC-T (ours)', dict(family='tecT', q=8.0, kappa=1.25, ramp=True, T0=20.0, S=S)),
    ('Onsager SCA (ours)', dict(family='onsager', q=8.0, lam=1.05, ramp=True, T0=12.0, S=S)),
]


def main():
    J, sumw = M.m.load()
    out = {}; arrays = {}
    for i, (name, cfg) in enumerate(METHODS):
        t = time.time(); E, cuts = M.run_method(J, sumw, cfg, 256, 71001 + i)
        arrays[name] = E.astype(np.float32)
        out[name] = dict(cfg=cfg, runs=256, mean_cut=float(cuts.mean()), sd_cut=float(cuts.std()), max_cut=float(cuts.max()),
                         p33000=float((cuts >= 33000).mean()), p33200=float((cuts >= 33200).mean()),
                         p33337=float((cuts >= 33337).mean()), final_E_mean=float(E[-1].mean()), sec=time.time() - t)
        print(name, json.dumps({k: v for k, v in out[name].items() if k != 'cfg'}), flush=True)
    np.savez_compressed(ROOT / 'trajectories_S960.npz', **{k.replace(' ', '_').replace('(', '').replace(')', ''): v for k, v in arrays.items()})
    (ROOT / 'final_S960.json').write_text(json.dumps(dict(S=S, sumw=sumw, methods=out), indent=1) + '\n')


if __name__ == '__main__':
    main()
