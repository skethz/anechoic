"""Precision-study analysis (PROTOCOL.md): per (problem, instance, rule, K) mean normalized quality, feasibility,
P(target), paired difference from full precision (95% normal interval over the 1024 paired trials), the smallest K
within tolerance (difference >= -0.005 at K and every larger tested K), and per-problem summaries.
Writes precision/precision_summary.json and precision/precision_analysis.txt."""
import json
import math
from collections import defaultdict
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
OUT = HERE / 'precision'
RULES = ('SCA', 'TEC', 'Onsager-kT', 'Onsager-online')
KS = (2, 3, 4, 6, 8)
TOL = 0.005


def main():
    runs = {}
    for f in (OUT / 'runs').glob('*.json'):
        d = json.loads(f.read_text())
        runs[(d['problem'], d['instance'], d['rule'], d['K'])] = d
    res = defaultdict(dict)
    for (prob, inst, rule, K), d in runs.items():
        if K is not None:
            continue
    keys = sorted({(p, i, r) for (p, i, r, K) in runs})
    summ = {}
    for p, i, r in keys:
        full = runs.get((p, i, r, None))
        if full is None:
            continue
        qf = np.array(full['quality'])
        row = dict(full=dict(mean_quality=float(qf.mean()), p_feasible=full['p_feasible'], p_target=full['p_target'],
                             S=full['S'], checks=full['checks']), K={})
        within = {}
        for K in KS:
            d = runs.get((p, i, r, K))
            if d is None:
                continue
            q = np.array(d['quality'])
            diff = q - qf
            se = float(diff.std(ddof=1) / math.sqrt(len(diff)))
            row['K'][K] = dict(mean_quality=float(q.mean()), p_feasible=d['p_feasible'], p_target=d['p_target'],
                               diff=float(diff.mean()), diff_ci95=(float(diff.mean() - 1.96 * se), float(diff.mean() + 1.96 * se)),
                               alpha=d['alpha'], checks=d['checks'], max_abs_b=d['qstats']['max_abs_b'],
                               bits_b=d['qstats']['bits_b_signed'], levels=d['qstats']['distinct_absJ_levels'])
            within[K] = diff.mean() >= -TOL
        kmin = None
        for K in KS:
            if all(within.get(k, False) for k in KS if k >= K):
                kmin = K; break
        row['K_min'] = kmin
        summ[f'{p}/{i}/{r}'] = row
    # per problem and rule
    agg = {}
    L = []
    for p in ('mcp', 'gpp', 'tsp'):
        L.append(f'\n## {p.upper()}: mean quality (feasible) by K, four engine rules; K_min = smallest K within -{TOL} of full\n')
        L.append('| Rule | full | ' + ' | '.join(f'K={K}' for K in KS) + ' | K_min max over instances | instances within tol. at K=2/3/4/6/8 |')
        L.append('|---|---|' + '---|' * len(KS) + '---|---|')
        for r in RULES:
            rows = [v for k, v in summ.items() if k.startswith(p + '/') and k.endswith('/' + r)]
            if not rows:
                continue
            fq = np.mean([v['full']['mean_quality'] for v in rows]); ff = np.mean([v['full']['p_feasible'] for v in rows])
            cells = []
            wcount = []
            for K in KS:
                ks = [v['K'][K] for v in rows if K in v['K']]
                cells.append(f"{np.mean([x['mean_quality'] for x in ks]):.4f} ({np.mean([x['p_feasible'] for x in ks]):.2f})")
                wcount.append(sum(1 for x in ks if x['diff'] >= -TOL))
            kmins = [v['K_min'] for v in rows]
            kmax = None if any(k is None for k in kmins) else max(kmins)
            agg[f'{p}/{r}'] = dict(n=len(rows), full_mean_quality=float(fq), full_p_feasible=float(ff),
                                   K_min_max=kmax, K_min_list=kmins, within_counts=dict(zip(KS, wcount)))
            L.append(f"| {r} | {fq:.4f} ({ff:.2f}) | " + ' | '.join(cells) +
                     f" | {kmax if kmax is not None else '>8 (some instance never within)'} | {'/'.join(map(str, wcount))} of {len(rows)} |")
    (OUT / 'precision_summary.json').write_text(json.dumps(dict(per_instance=summ, per_problem_rule=agg, tolerance=TOL),
                                                           indent=1, default=str) + '\n')
    (OUT / 'precision_analysis.txt').write_text('\n'.join(L) + '\n')
    print('\n'.join(L))


if __name__ == '__main__':
    main()
