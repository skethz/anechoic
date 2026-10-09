"""Where do the pilot selections of the four engine rules sit in their grids? (input to PROTOCOL amendment 1)
For each class, rule and S: share of selections on each grid edge, and the pilot landscape (best pilot mean cut per value
of each parameter, relative to the best overall, averaged over the class's instances). Pilot data only. edge_report.txt"""
import collections
import glob
import json
import sys
from pathlib import Path

sys.dont_write_bytecode = True
import run_gset as R  # noqa: E402

ROOT = Path(__file__).resolve().parent
SETS = {
    'SCA': lambda sh: dict(q=sh['plain_q'], T0=sh['plain_T0']),
    'TEC': lambda sh: dict(jv=sh['jv'], q=sh['q'], T0=sh['T0']),
    'Onsager-kT': lambda sh: dict(kappa=R.KAPPA, q=sh['q'], T0=sh['T0']),
    'Onsager-online': lambda sh: dict(lam_k2000=R.LAMBDA, q=sh['q'], T0=sh['T0']),
}


def main(only=None):
    L = []
    edge = collections.defaultdict(collections.Counter); n = collections.Counter()
    land = collections.defaultdict(lambda: collections.defaultdict(list))
    for f in sorted(glob.glob(str(ROOT / 'results/budget/G*/*.json'))):
        d = json.loads(Path(f).read_text())
        if d['method'] not in SETS or d['selected'] is None:
            continue
        g = int(d['instance'][1:])
        if only and g not in only:
            continue
        sets = SETS[d['method']](R.SHARED[d['wclass']])
        key = (d['wclass'], d['method'], d['S'])
        n[key] += 1
        rel = d['selected_cfg']['rel']
        for k, vals in sets.items():
            if rel[k] == min(vals):
                edge[key][f'{k}=min'] += 1
            if rel[k] == max(vals):
                edge[key][f'{k}=max'] += 1
        best = max(p['mean_cut'] for p in d['pilot'])
        for k, vals in sets.items():
            for v in vals:
                m = max(p['mean_cut'] for p in d['pilot'] if p['cfg']['rel'][k] == v)
                land[(d['wclass'], d['method'], k)][v].append(m - best)
    for key in sorted(n):
        L.append(f"{key[0]:4s} {key[1]:15s} S={key[2]:<5d} n={n[key]}: " + ', '.join(f'{e} {c}' for e, c in sorted(edge[key].items())))
    L.append('\nPilot landscape: mean over instances and S of (best pilot mean cut at this value - best overall), per parameter')
    for (c, m, k), d in sorted(land.items()):
        L.append(f"{c:4s} {m:15s} {k:10s}: " + '  '.join(f"{v:g}: {sum(x) / len(x):+.2f}" for v, x in sorted(d.items())))
    txt = '\n'.join(L)
    print(txt)
    return txt


if __name__ == '__main__':
    t = main()
    (ROOT / 'edge_report.txt').write_text(t + '\n')
