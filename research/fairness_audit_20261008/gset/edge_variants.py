"""Grid-edge counts of the G-set engine-rule selections, both grid variants (no new runs).
(a) per-budget pilot selections (51 instances x 5 budgets per rule); (b) the configuration measured on the board
(tts_E12 selection, one per instance and rule). An axis value counts as an edge if it equals the min or max of that
variant's value set for the rule and weight class (2-level axes are always at an edge). S* = 250 or 4000 is reported
separately (budget-set edge, identical for all rules)."""
import collections
import json
import sys
from pathlib import Path

GS = Path(__file__).resolve().parents[2] / 'gset_20261007'
sys.path.insert(0, str(GS))
sys.dont_write_bytecode = True
import run_gset as R  # noqa: E402
import run_gset_ext as X  # noqa: E402

PAR = {'SCA': None, 'TEC': 'jv', 'Onsager-kT': 'kappa', 'Onsager-online': 'lam_k2000'}


def value_sets(method, wc, variant):
    o = R.SHARED[wc]
    sets = {'SCA': dict(q=set(o['plain_q']), T0=set(o['plain_T0'])),
            'TEC': dict(jv=set(o['jv']), q=set(o['q']), T0=set(o['T0'])),
            'Onsager-kT': dict(kappa=set(R.KAPPA), q=set(o['q']), T0=set(o['T0'])),
            'Onsager-online': dict(lam_k2000=set(R.LAMBDA), q=set(o['q']), T0=set(o['T0']))}[method]
    if variant == 'extended_amendment1':
        e = X.EXT[wc]
        add = {'SCA': dict(q=e['plain_q'], T0=e['plain_T0']), 'TEC': dict(jv=e['jv'], q=e['q'], T0=e['T0']),
               'Onsager-kT': dict(kappa=X.KAPPA_EXT, q=e['q'], T0=e['T0']),
               'Onsager-online': dict(lam_k2000=X.LAMBDA_EXT, q=e['q'], T0=e['T0'])}[method]
        sets = {k: v | set(add[k]) for k, v in sets.items()}
    return sets


def edges(rel, sets):
    out = []
    for k, vals in sets.items():
        if rel[k] == min(vals): out.append(f'{k}=min')
        if rel[k] == max(vals): out.append(f'{k}=max')
    return out


def main():
    sc = json.loads((GS / 'selected_configs.json').read_text())['instances']
    res = {}
    lines = []
    for variant in ('original', 'extended_amendment1'):
        for scope in ('per_budget', 'tts_E12'):
            cnt = collections.defaultdict(collections.Counter); n = collections.Counter(); any_edge = collections.Counter()
            for gname, inst in sc.items():
                g = int(gname[1:]); wc = R.wclass(g)
                for rule in R.ENGINE_RULES:
                    v = inst['rules'][rule][variant]
                    cfgs = list(v['per_budget'].values()) if scope == 'per_budget' else [v['tts_E12']['cfg']]
                    sets = value_sets(rule, wc, variant)
                    for cfg in cfgs:
                        e = edges(cfg['rel'], sets)
                        for k, vals in sets.items():
                            cnt[(rule, 'all')][f'{k}:levels={len(vals)}'] += 0
                        for x in e:
                            cnt[(rule, 'all')][x] += 1; cnt[(rule, wc)][x] += 1
                        if scope == 'tts_E12':
                            if cfg['S'] == 250: cnt[(rule, 'all')]['S*=250(min budget)'] += 1
                            if cfg['S'] == 4000: cnt[(rule, 'all')]['S*=4000(max budget)'] += 1
                        any_edge[(rule, 'all')] += bool([x for x in e if not x.startswith('q') or True])
                        n[(rule, 'all')] += 1; n[(rule, wc)] += 1
            key = f'{variant}/{scope}'
            res[key] = {f'{r}/{w}': dict(n=n[(r, w)], edges=dict(cnt[(r, w)])) for (r, w) in sorted(n)}
            lines.append(f'=== {key}')
            for (r, w) in sorted(n):
                c = cnt[(r, w)]
                lines.append(f"{r:15s} {w:4s} n={n[(r, w)]:3d}: " + ', '.join(
                    f"{k} {v} ({100 * v / n[(r, w)]:.0f}%)" for k, v in sorted(c.items()) if v))
    t = '\n'.join(lines)
    print(t)
    Path(__file__).with_name('edge_variants.txt').write_text(t + '\n')
    Path(__file__).with_name('edge_variants.json').write_text(json.dumps(res, indent=1) + '\n')


if __name__ == '__main__':
    main()
