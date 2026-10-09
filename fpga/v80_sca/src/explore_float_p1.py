"""Exploratory (not part of PROTOCOL_HW): fresh float-model cohorts for P1/P2/T1 with the J sweep's own kernel (joint.run_sca)."""
import json, sys
from pathlib import Path
import numpy as np
sys.path.insert(0, str(__import__('pathlib').Path(__file__).resolve().parent / '../../../research/theory_ideas_20261003'))
import sweeps_abc as m
import joint
J, sumw = m.load()[:2]
B = int(sys.argv[1]) if len(sys.argv) > 1 else 4096
cfgs = [('P1', 'plain', 4.0, 0.0, 30.0, 1560), ('P2', 'plain', 6.0, 0.0, 30.0, 1560), ('T1', 'tec', 8.0, -4.0, 30.0, 1560)]
out = {}
for (key, mode, q, par, t0, S), ss in zip(cfgs, np.random.SeedSequence(20261004_7).spawn(len(cfgs))):
    rng = np.random.default_rng(ss); k = 0; n = 0
    for _ in range(B // 512):
        s0 = rng.choice(np.array([-1.0, 1.0], dtype=np.float32), size=(512, len(J)))
        s, fl = joint.run_sca(J, s0, S, t0, q, mode, par, False, rng)
        c = m.cut(J, sumw, s); k += int((c >= m.TARGET).sum()); n += 512
    out[key] = dict(trials=n, successes=k, p=k / n)
    print(key, json.dumps(out[key]), flush=True)
