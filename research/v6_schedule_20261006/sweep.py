"""Schedule re-selection for the v6 engine with 9 concurrent engines (exploratory pilot -> fresh held-out check).

Why: the frozen O1-O4 were chosen with the v5 cost model (31 cycles/step) for E = 1/4/16. v6 costs
0.2714*flips + 16.91*S + 1148 cycles per trial (RTL simulation fit) plus S + 300 cycles of per-launch streaming, and the
9-engine primary estimator (Eq. 9 on measured rounds, P_round = 1 -> one round) rewards the fastest schedule whose round
success is about 0.99. Same CPU model as the ablations (research/ablation_20261005/abl.py: clipped parallel SCA, float32).

Stages: pilot (256 runs/config, seed 61001) -> selection by Wilson lower bound -> held-out (2048 fresh runs, seed 61002).
"""
import itertools
import json
import math
import multiprocessing as mp
import sys
import time
from pathlib import Path

import numpy as np

sys.path.insert(0, str(__import__('pathlib').Path(__file__).resolve().parent / '../ablation_20261005'))
import abl  # noqa: E402

ROOT = Path(__file__).resolve().parent
E = 9
MHZ = 300.0
ROUND_MAX = 1.035   # mean of max over 9 engines / mean launch (from the known-trial prediction, O4/O3)


def t_round_ms(flips, S):
    return ROUND_MAX * (0.2714 * flips + 16.91 * S + 1148 + S + 300) / (MHZ * 1e3)


def r99(p):
    return abl.r99(p)


def tts_primary(t, p):
    P = 1 - (1 - p) ** E
    return t if P >= 0.995 else t * r99(P)


def tts_secondary(t, p):
    if p <= 0: return math.inf
    if p >= 1: return t
    return t * math.log(0.01) / (E * math.log1p(-p))


def grid():
    g = []
    for q, lam, T0, S in itertools.product((4.0, 6.0, 8.0), (0.8, 0.9, 1.05), (12.0, 15.0, 20.0), (200, 240, 280, 320, 360)):
        g.append(dict(family='onsager', q=q, lam=lam, ramp=True, T0=T0, S=S))
    for q, kap, T0, S in itertools.product((6.0, 8.0), (1.0, 1.25, 1.5, 1.75), (15.0, 20.0, 25.0), (200, 240, 280, 320, 360)):
        g.append(dict(family='tecT', q=q, kappa=kap, ramp=True, T0=T0, S=S))
    for q, kap, T0, S in itertools.product((6.0, 8.0), (1.25, 1.5, 1.75, 2.0), (15.0, 20.0, 25.0), (560, 760, 960, 1260)):
        g.append(dict(family='tecT', q=q, kappa=kap, ramp=True, T0=T0, S=S))
    for q, lam, T0, S in itertools.product((6.0, 8.0), (0.9, 1.05, 1.2), (12.0, 15.0, 20.0), (560, 760, 960, 1260)):
        g.append(dict(family='onsager', q=q, lam=lam, ramp=True, T0=T0, S=S))
    return g


_J = None


def _init():
    global _J
    _J = abl.m.load()


def _eval(args):
    cfg, B, seed = args
    J, sumw = _J
    rng = np.random.default_rng(seed); t1 = time.perf_counter()
    s0 = rng.choice(np.array([-1.0, 1.0], dtype=np.float32), size=(B, len(J)))
    s, fl = abl.run(J, s0, cfg, rng)
    c = abl.m.cut(J, sumw, s); k = int((c >= abl.m.TARGET).sum())
    t = t_round_ms(float(fl.mean()), cfg['S'])
    return dict(cfg, trials=B, successes=k, p=k / B, mean_flips=float(fl.mean()), t_round_ms=t,
                tts_primary_ms=tts_primary(t, k / B), tts_secondary_ms=tts_secondary(t, k / B), sec=time.perf_counter() - t1)


def run_all(cfgs, B, seed_root, tag):
    seeds = np.random.SeedSequence(seed_root).spawn(len(cfgs))
    with mp.Pool(processes=max(1, mp.cpu_count() - 2), initializer=_init) as pool:
        out = []
        for r in pool.imap(_eval, [(c, B, s) for c, s in zip(cfgs, seeds)]):
            out.append(r)
            print(tag, json.dumps({k: (round(v, 4) if isinstance(v, float) else v) for k, v in r.items() if k != 'sec'}), flush=True)
    return out


def select(rows, n_per=3):
    """Per family and objective: rank by the TTS at the Wilson lower bound of p (conservative), keep the top n."""
    sel = []
    for fam in ('onsager', 'tecT'):
        for obj, f in (('primary', tts_primary), ('secondary', tts_secondary)):
            cand = [r for r in rows if r['family'] == fam]
            cand.sort(key=lambda r: f(r['t_round_ms'], abl.wilson_lower(r['successes'], r['trials'])))
            for r in cand[:n_per]:
                c = {k: r[k] for k in ('family', 'q', 'lam', 'kappa', 'ramp', 'T0', 'S') if k in r}
                if c not in [s['cfg'] for s in sel]:
                    sel.append(dict(cfg=c, objective=obj))
    return sel


def main():
    stage = sys.argv[1] if len(sys.argv) > 1 else 'all'
    if stage in ('pilot', 'all'):
        pilot = run_all(grid(), 256, 61001, 'PILOT')
        (ROOT / 'pilot.json').write_text(json.dumps(pilot, indent=1) + '\n')
    if stage in ('holdout', 'all'):
        pilot = json.loads((ROOT / 'pilot.json').read_text())
        sel = select(pilot)
        hold = run_all([s['cfg'] for s in sel], 2048, 61002, 'HOLD')
        for s, h in zip(sel, hold): h['objective'] = s['objective']
        (ROOT / 'holdout.json').write_text(json.dumps(dict(selection=sel, holdout=hold), indent=1) + '\n')


if __name__ == '__main__':
    main()
