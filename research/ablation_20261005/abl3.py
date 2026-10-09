"""ABL3 (PROTOCOL_ABL3.md): short schedules for TEC-T and Onsager; selection over all pilots; paired held-out (seed 50002)."""
import itertools
import json
from pathlib import Path

import numpy as np

import abl

ROOT = Path(__file__).resolve().parent
J, sumw = abl.m.load()
cfgs = [dict(family='tecT', q=q, kappa=k, ramp=r, T0=T0, S=S)
        for S, k, r, q, T0 in itertools.product((200, 280, 360), (0.75, 1.0, 1.25, 1.5), (False, True), (6.0, 8.0), (12.0, 15.0, 20.0))]
cfgs += [dict(family='onsager', q=q, lam=l, ramp=r, T0=T0, S=S)
         for S, l, r, q, T0 in itertools.product((200, 280), (0.9, 1.05), (False, True), (6.0, 8.0), (12.0, 15.0))]
pilot = abl.evaluate(J, sumw, cfgs, 256, np.random.SeedSequence(50021).spawn(len(cfgs)), tag='PILOT3')
(ROOT / 'pilot_abl3.json').write_text(json.dumps(pilot, indent=1) + '\n')
tecT = [r for r in json.loads((ROOT / 'pilot_new.json').read_text()) if r['family'] == 'tecT']
tecT += json.loads((ROOT / 'pilot_tecT2.json').read_text()) + [r for r in pilot if r['family'] == 'tecT']
ons = [r for r in abl.j_pilot_as_cfgs() if r['family'] == 'onsager'] + [r for r in pilot if r['family'] == 'onsager']
sel = abl.select(tecT + ons)
prior = json.loads((ROOT / 'holdout.json').read_text())['holdout'] + json.loads((ROOT / 'holdout_tecT2.json').read_text())['holdout']
KEYS = ('family', 'q', 'T0', 'S', 'lam', 'ramp', 'kappa')
done = {tuple((k, r.get(k)) for k in KEYS) for r in prior if r['family'] in ('tecT', 'onsager')}
todo = []
for v in sel.values():
    kk = tuple((k, v.get(k)) for k in KEYS)
    if kk not in done and v not in todo: todo.append(v)
hold = abl.evaluate(J, sumw, todo, 1024, [np.random.SeedSequence(50002)] * len(todo), tag='HOLD3') if todo else []
(ROOT / 'holdout_abl3.json').write_text(json.dumps(dict(selection=sel, holdout=hold), indent=1) + '\n')
