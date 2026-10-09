"""ABL2 (PROTOCOL_ABL2.md): extended TEC-T grid; selection over TEC-T + TEC-T2 pilots; paired held-out (seed 50002)."""
import itertools
import json
import sys
from pathlib import Path

import numpy as np

import abl

ROOT = Path(__file__).resolve().parent
J, sumw = abl.m.load()
stage = sys.argv[1] if len(sys.argv) > 1 else 'all'
if stage in ('pilot', 'all'):
    cfgs = [dict(family='tecT', q=q, kappa=k, ramp=r, T0=T0, S=S)
            for k, r, q, T0, S in itertools.product((0.75, 1.0, 1.25, 1.5, 2.0), (False, True), (4.0, 6.0, 8.0), (20.0, 25.0, 30.0), (560, 960, 1560))]
    pilot = abl.evaluate(J, sumw, cfgs, 256, np.random.SeedSequence(50011).spawn(len(cfgs)), tag='PILOT2')
    (ROOT / 'pilot_tecT2.json').write_text(json.dumps(pilot, indent=1) + '\n')
if stage in ('holdout', 'all'):
    rows = [r for r in json.loads((ROOT / 'pilot_new.json').read_text()) if r['family'] == 'tecT']
    rows += json.loads((ROOT / 'pilot_tecT2.json').read_text())
    sel = abl.select(rows)
    distinct = []
    for v in sel.values():
        if v not in distinct: distinct.append(v)
    hold = abl.evaluate(J, sumw, distinct, 1024, [np.random.SeedSequence(50002)] * len(distinct), tag='HOLD2')
    (ROOT / 'holdout_tecT2.json').write_text(json.dumps(dict(selection=sel, holdout=hold), indent=1) + '\n')
