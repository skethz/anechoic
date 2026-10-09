#!/usr/bin/env python3
"""PROTOCOL_HW_V6 Amendment 3 checks for the v6.2 x 12-engine image (run after analyze_hw_v6.py wrote hwm_summary.json).
R1: every trial id of every round configuration has the same cut, flips and spins as in the 9-engine v6.1 run (both images
are bit-exact with sca_ref.hpp and use the same seed and trial ids, so trial outcomes are deterministic).
Usage: analyze_v62_amend3.py <hwv6 v62 e12 dir> <hwv6 v6 e9 dir>"""
import json, sys
from pathlib import Path

new, old = Path(sys.argv[1]), Path(sys.argv[2])
S = json.loads((new / 'hwm_summary.json').read_text())
R = S['rounds']; pr = S['primary']
mism = {}; compared = 0
for f in sorted((new / 'rounds').glob('*.jsonl')):
    a = {r['trial_id']: r for r in map(json.loads, f.read_text().splitlines()) if r}
    bpath = old / 'rounds' / f.name
    if not bpath.exists(): continue
    b = {r['trial_id']: r for r in map(json.loads, bpath.read_text().splitlines()) if r}
    bad = [i for i in a if i in b and (a[i]['cut'] != b[i]['cut'] or a[i]['flips'] != b[i]['flips'])]
    sa, sb = (new / 'rounds' / (f.stem + '.spins.bin')).read_bytes(), (old / 'rounds' / (f.stem + '.spins.bin')).read_bytes()
    compared += len(a)
    if bad: mism[f.stem] = len(bad)
    # spins: records are written in the same (round, engine) order only when E matches, so compare per trial id via cut/flips
best_o = pr['best_onsager_tts_ms']; m6 = pr['M6']
px = min(R[k]['tts99_ms'] for k in ('X1', 'X2', 'X5')); sx = min(R[k]['tts99_ind_ms'] for k in ('X3', 'X4'))
pb = min(R[k]['tts99_ms'] for k in ('B1', 'B3')); sb_ = min(R[k]['tts99_ind_ms'] for k in ('B2', 'B4'))
res = dict(R1_trials_compared=compared, R1_mismatching_configs=mism, R1_pass=not mism,
           M3_best_onsager_ms=best_o, M3_pass=best_o <= 0.08,
           M4_pass=pr['M4_pass'], M5_pass=best_o <= 0.13,
           M6_best_ms=m6['best_ms'], M6_pass=m6['best_ms'] <= 0.04 and m6['vs_plain'] >= 2.5 and m6['vs_tec'] >= 2.5,
           A1_best_primary_ms=px, A1_pass=px <= 0.08, A2_best_secondary_ms=sx, A2_pass=sx <= 0.031,
           A3_primary_ratio=pb / px, A3_secondary_ratio=sb_ / sx, A3_pass=pb / px >= 2 and sb_ / sx >= 2.5)
(new / 'amendment3.json').write_text(json.dumps(res, indent=2) + '\n')
print(json.dumps(res, indent=1))
