"""Remaining grid edges after PROTOCOL_AMENDMENT1.md: share of union selections at the extreme values of the union value
sets (per class and engine rule, all instances and budgets). Pilot-based selections only. Output: edge_report_ext.txt"""
import collections
import glob
import json
import sys
from pathlib import Path

sys.dont_write_bytecode = True
import run_gset as R  # noqa: E402
import run_gset_ext as X  # noqa: E402

ROOT = Path(__file__).resolve().parent
PAR = {'SCA': None, 'TEC': 'jv', 'Onsager-kT': 'kappa', 'Onsager-online': 'lam_k2000'}


def value_sets(method, wc):
    o, e = R.SHARED[wc], X.EXT[wc]
    if method == 'SCA':
        return dict(q=set(o['plain_q']) | set(e['plain_q']), T0=set(o['plain_T0']) | set(e['plain_T0']))
    p = {'TEC': (o['jv'], e['jv']), 'Onsager-kT': (R.KAPPA, X.KAPPA_EXT), 'Onsager-online': (R.LAMBDA, X.LAMBDA_EXT)}[method]
    return {PAR[method]: set(p[0]) | set(p[1]), 'q': set(o['q']) | set(e['q']), 'T0': set(o['T0']) | set(e['T0'])}


def main():
    cnt = collections.defaultdict(collections.Counter); n = collections.Counter(); fromext = collections.Counter()
    for f in glob.glob(str(ROOT / 'results/ext_budget/G*/*.json')):
        d = json.loads(Path(f).read_text())
        k = (d['wclass'], d['method']); n[k] += 1; fromext[k] += d['from_extension']
        rel = d['union_selected_cfg']['rel']
        for par, vals in value_sets(d['method'], d['wclass']).items():
            if rel[par] == min(vals):
                cnt[k][f'{par}=min({min(vals):g})'] += 1
            if rel[par] == max(vals):
                cnt[k][f'{par}=max({max(vals):g})'] += 1
    L = ['Union-grid selections (all instances x 5 budgets) at the extremes of the union value sets']
    for k in sorted(n):
        L.append(f"{k[0]:4s} {k[1]:15s} n={n[k]:3d} from extension {fromext[k]:3d}: " +
                 ', '.join(f'{e} {c}' for e, c in sorted(cnt[k].items())))
    t = '\n'.join(L)
    print(t)
    (ROOT / 'edge_report_ext.txt').write_text(t + '\n')


if __name__ == '__main__':
    main()
