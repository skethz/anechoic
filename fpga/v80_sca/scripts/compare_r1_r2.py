#!/usr/bin/env python3
"""PROTOCOL_HW_MB Amendment 2: compare the revision-2 image's runs with the first multi-bit image's, trial id by trial id.
Same algorithm and arithmetic (only the device's cut word changed), same seed, trial ids and launch shape, so every trial
must have identical spins, flips and independently rescored cut; on r2, the device cut must equal the rescored cut for
every trial (on r1 it did not when sum s*h < 0).
Usage: compare_r1_r2.py <r1 result dir> <r2 result dir> [subdirs ...]   (default subdirs: rounds gset)
Writes <r2 dir>/r1_r2_compare.json."""
import json
import sys
from pathlib import Path


def load(p):
    return [json.loads(line) for line in Path(p).read_text().splitlines() if line.strip()]


def main():
    r1, r2 = Path(sys.argv[1]), Path(sys.argv[2]); subs = sys.argv[3:] or ['rounds', 'gset']
    out = {}; total = dict(trials=0, mismatched=0, r2_disagree=0, r1_disagree=0, files=0, missing=0)
    for sub in subs:
        for f2 in sorted((r2 / sub).glob('*.jsonl')):
            f1 = r1 / sub / f2.name
            if not f1.exists(): total['missing'] += 1; continue
            a, b = load(f1), load(f2)
            s1, s2 = (r1 / sub / (f2.stem + '.spins.bin')).read_bytes(), (r2 / sub / (f2.stem + '.spins.bin')).read_bytes()
            ia = {r['trial_id']: (i, r) for i, r in enumerate(a)}
            mism = 0
            for j, r in enumerate(b):
                i, q = ia.get(r['trial_id'], (None, None))
                if q is None or q['flips'] != r['flips'] or q['cut'] != r['cut'] or s1[256 * i:256 * i + 256] != s2[256 * j:256 * j + 256]:
                    mism += 1
            d2 = sum(not r['scores_agree'] for r in b); d1 = sum(not r['scores_agree'] for r in a)
            out[f'{sub}/{f2.stem}'] = dict(trials=len(b), mismatched=mism, r2_scores_disagree=d2, r1_scores_disagree=d1)
            total['trials'] += len(b); total['mismatched'] += mism; total['r2_disagree'] += d2; total['r1_disagree'] += d1; total['files'] += 1
    res = dict(total=total, per_file=out, A2_R_pass=total['files'] > 0 and total['mismatched'] == 0 and total['missing'] == 0,
               A2_C_pass=total['files'] > 0 and total['r2_disagree'] == 0)
    (r2 / 'r1_r2_compare.json').write_text(json.dumps(res, indent=1) + '\n')
    print(json.dumps(total), 'A2-R', res['A2_R_pass'], 'A2-C', res['A2_C_pass'])
    for k, v in out.items():
        if v['mismatched'] or v['r2_scores_disagree'] or v['r1_scores_disagree']: print(k, v)


if __name__ == '__main__':
    main()
