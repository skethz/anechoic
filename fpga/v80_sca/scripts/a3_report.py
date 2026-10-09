#!/usr/bin/env python3
"""Paper numbers from an Amendment-3 summary (a3_summary.json); reporting only, no verdicts.
- Geometric-mean primary TTS99 per rule (SCA, TEC, Onsager-kT, Onsager-online) over all instances and per graph class
  (random +1, random +-1, toroidal +-1, planar +1, planar +-1), separately for the pre-registered ('original') grid and the
  'extended_amendment1' grid (a cohort shared by both variants counts in both).
- Instances where Onsager-kT's TTS99 is below plain SCA's, below TEC's, and below both (original grid).
Weight signs come from the G-set files themselves (any negative weight -> +-1 class).
Usage: a3_report.py <a3_summary.json> <G-set data dir> [--md]"""
import json
import math
import sys
from pathlib import Path

RULES = ['SCA', 'TEC', 'Onsager-kT', 'Onsager-online']
GRIDS = [('original', 'pre-registered (original) grid'), ('extended_amendment1', 'extended grid (extended_amendment1)')]
CLASSES = ['random +1', 'random ±1', 'toroidal ±1', 'planar +1', 'planar ±1']


def topology(i):
    if i <= 10 or 22 <= i <= 31 or 43 <= i <= 47: return 'random'
    if 11 <= i <= 13 or 32 <= i <= 34: return 'toroidal'
    return 'planar'


def gmean(xs):
    xs = list(xs)
    if not xs: return None, 0, 0
    inf = sum(1 for x in xs if not math.isfinite(x))
    fin = [x for x in xs if math.isfinite(x) and x > 0]
    if inf: return math.inf, len(xs), inf
    return math.exp(sum(map(math.log, fin)) / len(fin)), len(xs), 0


def main():
    S = json.load(open(sys.argv[1])); gdir = Path(sys.argv[2]); md = '--md' in sys.argv
    C = {k: c for k, c in S['cohorts'].items() if isinstance(c, dict)}
    insts = sorted({c['instance'] for c in C.values()}, key=lambda s: int(s[1:]))
    cls = {}
    for i in insts:
        neg = False
        with open(gdir / i) as f:
            next(f)
            for line in f:
                p = line.split()
                if len(p) >= 3 and int(p[2]) < 0: neg = True; break
        cls[i] = f"{topology(int(i[1:]))} {'±1' if neg else '+1'}"
    tts = {}   # (grid, rule, instance) -> tts99
    for c in C.values():
        for g, _ in GRIDS:
            if g in c['variants']:
                tts[(g, c['rule'], c['instance'])] = c['tts99_ms']
    out = dict(classes={i: cls[i] for i in insts}, geomeans={}, onsager_kT_wins={})
    for g, gname in GRIDS:
        rows = []
        for r in RULES:
            row = {}
            for cl in CLASSES + ['all']:
                xs = [tts[(g, r, i)] for i in insts if (cl == 'all' or cls[i] == cl) and (g, r, i) in tts]
                gm, n, ninf = gmean(xs)
                row[cl] = dict(geomean_tts99_ms=gm, instances=n, no_success=ninf)
            out['geomeans'][f'{g}/{r}'] = row
            rows.append((r, row))
        if md:
            print(f"\n**Geometric-mean primary TTS99 (ms), {gname}, 12 engines.** In brackets: the number of instances.\n")
            print('| Rule | ' + ' | '.join(CLASSES) + ' | all |')
            print('|---|' + '---:|' * (len(CLASSES) + 1))
            for r, row in rows:
                cells = []
                for cl in CLASSES + ['all']:
                    v = row[cl]
                    cells.append('–' if v['geomean_tts99_ms'] is None else
                                 (f"∞ ({v['instances']}, {v['no_success']} without success)" if not math.isfinite(v['geomean_tts99_ms'])
                                  else f"{v['geomean_tts99_ms']:.3f} ({v['instances']})"))
                print(f'| {r} | ' + ' | '.join(cells) + ' |')
    g = 'original'
    w = dict(vs_SCA=[], vs_TEC=[], vs_both=[], compared=0)
    for i in insts:
        k, s, t = tts.get((g, 'Onsager-kT', i)), tts.get((g, 'SCA', i)), tts.get((g, 'TEC', i))
        if k is None or s is None or t is None: continue
        w['compared'] += 1
        if k < s: w['vs_SCA'].append(i)
        if k < t: w['vs_TEC'].append(i)
        if k < s and k < t: w['vs_both'].append(i)
    out['onsager_kT_wins'] = w
    print(f"\nOnsager-kT TTS99 below plain SCA on {len(w['vs_SCA'])} of {w['compared']} instances, below TEC on {len(w['vs_TEC'])}, "
          f"below both on {len(w['vs_both'])} (original grid).")
    print('  not below SCA:', [i for i in insts if i not in w['vs_SCA']])
    print('  not below TEC:', [i for i in insts if i not in w['vs_TEC']])
    Path(sys.argv[1]).with_name('a3_report.json').write_text(json.dumps(out, indent=1) + '\n')


if __name__ == '__main__':
    main()
