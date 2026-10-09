"""Analysis of PROTOCOL_PUB.md (run after run_pub.py summarize): claim checks, tables, and the Figure 2 data in the exact
format and method names of research/algorithm_compare_20261007/plot_alg.py (written to fig/results/S<S>.json, so that an
unmodified copy of plot_alg.py in fig/ draws the new figure). Independent recomputation of p, MCS99 and mean energy from
the per-run final cuts stored in the results (not from the stored summary fields)."""
import json
import math
import shutil
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
RES = HERE / 'results_pub'
FIG = HERE / 'fig'
ROOT = HERE.parents[2]
S_LIST = (250, 500, 1000, 2000, 4000)
SUMW = -1040.0
# new name -> Figure 2 name (plot_alg.py ORDER)
FIGNAME = {'SA (Neal)': 'SA', 'SCA (STATICA)': 'SCA', 'TEC (published)': 'TEC', 'APC-SCA (published)': 'APC-SCA',
           'ReAIM ASA': 'ReAIM ASA', 'aSB': 'aSB', 'bSB': 'bSB', 'dSB': 'dSB', 'TEC-T (ours)': 'TEC-T (ours)',
           'Onsager SCA (ours)': 'Onsager SCA (ours)'}
SB = ('aSB', 'bSB', 'dSB')
OURS = ('TEC-T (ours)', 'Onsager SCA (ours)')


def mcs99(S, p):
    return math.inf if p <= 0 else (S if p >= 1 else S * math.log(0.01) / math.log1p(-p))


def recompute(f, S):
    c = np.array(f['cuts'], float)
    p = float((c >= 33000).mean())
    return dict(mean_cut=float(c.mean()), mean_H=float(SUMW - 2 * c.mean()), p=p, mcs99=mcs99(S, p), n=len(c))


def main():
    res = {S: json.loads((RES / f'S{S}.json').read_text()) for S in S_LIST}
    summ = json.loads((RES / 'summary.json').read_text())
    L = []
    rec = {S: {n: recompute(r['final'], S) for n, r in res[S].items()} for S in S_LIST}
    for S in S_LIST:   # stored vs recomputed
        for n, r in res[S].items():
            f = r['final']; x = rec[S][n]
            assert (f['mean_cut'] == x['mean_cut'] or abs(f['mean_cut'] - x['mean_cut']) < 1e-9) and abs(f['p33000'] - x['p']) < 1e-12, (S, n)
    L.append('Mean final cut / P(cut >= 33,000) / MCS99, final state, 256 runs (ours: frozen finals)')
    names = list(res[S_LIST[0]].keys())
    for S in S_LIST:
        L.append(f'--- S = {S}')
        order = sorted(names, key=lambda n: -rec[S][n]['mean_cut'])
        for n in order:
            x = rec[S][n]; bv = res[S][n]['final'].get('best_visited', {})
            L.append(f"{n:22s} mean cut {x['mean_cut']:9.1f}  p {x['p']:.3f}  MCS99 {x['mcs99']:9.0f}"
                     + (f"  | best-visited p {bv['p33000']:.3f}" if bv else '') + f"  cfg {res[S][n]['final'].get('cfg')}")
    # claim (a): lowest mean energy among discrete methods at every budget
    L.append('\nClaim (a): among discrete methods (all but aSB, bSB, dSB), do both forms reach lower mean energy than every baseline?')
    ok_all = True
    for S in S_LIST:
        disc = [n for n in names if n not in SB]
        best_base = max((n for n in disc if n not in OURS), key=lambda n: rec[S][n]['mean_cut'])
        both = all(rec[S][o]['mean_cut'] > rec[S][best_base]['mean_cut'] for o in OURS)
        ok_all &= both
        L.append(f"S={S}: best discrete baseline {best_base} {rec[S][best_base]['mean_cut']:.1f}; ours "
                 f"{rec[S]['TEC-T (ours)']['mean_cut']:.1f} / {rec[S]['Onsager SCA (ours)']['mean_cut']:.1f} -> {'YES' if both else 'NO'}"
                 f"; overall best {max(names, key=lambda n: rec[S][n]['mean_cut'])}")
    L.append(f'Claim (a) holds at every budget: {ok_all}')
    # claim (b): S = 500 success
    L.append('\nClaim (b): P(cut >= 33,000) at S = 500: ' + ', '.join(f"{n} {rec[500][n]['p']:.3f}" for n in names))
    # claim (c): steps to solution
    L.append('\nClaim (c): steps to solution, min over S (final state)')
    mins = {n: min((rec[S][n]['mcs99'], S) for S in S_LIST) for n in names}
    for n in sorted(names, key=lambda n: mins[n][0]):
        L.append(f'{n:22s} {mins[n][0]:9.0f} at S = {mins[n][1]}')
    disc_rank = sorted((n for n in names if n not in SB), key=lambda n: mins[n][0])
    L.append(f'Most step-efficient discrete method: {disc_rank[0]}; overall: {min(names, key=lambda n: mins[n][0])}')
    # sensitivity: the non-selected published values
    L.append('\nSensitivity: all published values (finals), min over S of MCS99 and mean cut at S = 1000')
    for n, per in summ['methods'].items():
        for k in sorted({int(k) for S in per for k in per[str(S)]['configs']} if per else []):
            vals = []
            for S in S_LIST:
                c = per.get(str(S), {}).get('configs', {}).get(str(k))
                if c:
                    vals.append((c['mcs99'], S, c['mean_cut'], c['p33000']))
            if vals:
                L.append(f"{n:22s} k={k} cfg@1000={per['1000']['configs'][str(k)]['cfg'] if '1000' in per else ''}: "
                         f"min MCS99 {min(vals)[0]:.0f} (S={min(vals)[1]}); p by S " + ' '.join(f'{v[3]:.2f}' for v in vals)
                         + f"; selected k={[per[str(S)]['selected_k'] for S in S_LIST if str(S) in per]}")
    t = '\n'.join(L)
    print(t)
    (RES / 'claims.txt').write_text(t + '\n')
    # Figure 2 data for an unmodified copy of plot_alg.py
    (FIG / 'results').mkdir(parents=True, exist_ok=True)
    shutil.copy(ROOT / 'research/algorithm_compare_20261007/plot_alg.py', FIG / 'plot_alg.py')
    for S in S_LIST:
        out = {FIGNAME[n]: dict(final=res[S][n]['final']) for n in res[S] if n in FIGNAME}
        (FIG / 'results' / f'S{S}.json').write_text(json.dumps(out) + '\n')
    print('figure data written to', FIG / 'results')


if __name__ == '__main__':
    main()
