"""Summaries of the exploratory calibration (calib.py). Quality per run: ref/value if feasible else 0, where ref is the
best feasible value found on that synthetic instance by any calibration point (GPP) or the exact optimum (TSP).
Usage: python3 calib_analyze.py stage1|stage2|penalty [top]"""
import json
import sys
from collections import defaultdict
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
OUT = HERE / 'calib'


def load(stage):
    recs = defaultdict(list)
    for f in sorted((OUT / stage).rglob('*.json')):
        d = json.loads(f.read_text())
        recs[(d['problem'], d['instance'])].append(d)
    return recs


def quality(prob, pt, ref):
    if prob == 'gpp':
        cut = np.array(pt['cut']); M = np.array(pt['M']); ok = np.array(pt['ok'])
        feas = (M == 0) & ok
        q = np.where(feas, ref / np.maximum(cut, 1), 0.0)
        return q, feas
    L = np.array(pt['L']); feas = np.array(pt['feasible'])
    return np.where(feas, ref / np.maximum(L, 1), 0.0), feas


def summarize(stage, top=5, refs=None):
    recs = load(stage)
    refs = dict(refs or {})
    best = {}
    for (prob, inst), ds in recs.items():
        if prob == 'gpp' and inst not in refs:
            m = np.inf
            for d in ds:
                for pt in d['points']:
                    cut = np.array(pt['cut']); feas = (np.array(pt['M']) == 0) & np.array(pt['ok'])
                    if feas.any():
                        m = min(m, cut[feas].min())
            refs[inst] = float(m)
        elif prob == 'tsp':
            refs[inst] = float(ds[0]['ref'])
        by_m = defaultdict(list)
        for d in ds:
            by_m[d['method'].split('#')[0]].extend(d['points'])
        print(f'\n=== {stage} {inst} ref={refs[inst]}')
        for meth, pts in sorted(by_m.items()):
            rows = []
            for pt in pts:
                q, feas = quality(prob, pt, refs[inst])
                rows.append((float(q.mean()), float(feas.mean()), pt['cfg']))
            rows.sort(key=lambda r: -r[0])
            best[(inst, meth)] = rows[0]
            print(f'  {meth:15s} n={len(rows):3d} best q={rows[0][0]:.4f} feas={rows[0][1]:.2f}')
            for r in rows[:top]:
                c = r[2]
                rel = c.get('rel', {k: c[k] for k in ('dt', 'xi') if k in c})
                print(f'      q={r[0]:.4f} feas={r[1]:.2f} {json.dumps(rel)}')
    return best, refs


if __name__ == '__main__':
    st = sys.argv[1]; top = int(sys.argv[2]) if len(sys.argv) > 2 else 5
    refs = None
    if st != 'stage1' and (OUT / 'refs_stage1.json').exists():
        refs = json.loads((OUT / 'refs_stage1.json').read_text())
    best, refs = summarize(st, top, refs)
    if st == 'stage1':
        (OUT / 'refs_stage1.json').write_text(json.dumps(refs, indent=1) + '\n')
