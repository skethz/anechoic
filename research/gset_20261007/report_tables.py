"""Markdown tables for REPORT.md from results_summary[_ext].json (no numbers are typed by hand).
Usage: python3 report_tables.py [results_summary.json | results_summary_ext.json]"""
import json
import math
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent
D = json.loads((ROOT / (sys.argv[1] if len(sys.argv) > 1 else 'results_summary.json')).read_text())
INST = D['instances']
RULES = ('SCA', 'TEC', 'Onsager-kT', 'Onsager-online')
SH = {'SCA': 'plain', 'TEC': 'TEC', 'Onsager-kT': 'kT', 'Onsager-online': 'online'}
METHODS = ('SA', 'SCA', 'TEC', 'APC-SCA', 'ReAIM ASA', 'aSB', 'bSB', 'dSB', 'Onsager-kT', 'Onsager-online')
CLASS_ORDER = ['R800+', 'R800+-', 'T800+-', 'P800+', 'P800+-', 'R2000+', 'R2000+-', 'T2000+-', 'P2000+', 'P2000+-', 'R1000+', 'P1000+']


def f(x, d=3):
    if x is None:
        return '–'
    if isinstance(x, (int,)) and not isinstance(x, bool):
        return str(x)
    if math.isinf(x):
        return '∞'
    if x == 0:
        return '0'
    return f'{x:.{d}g}'


def ms(x):
    if x is None:
        return '–'
    if math.isinf(x):
        return '∞'
    return f'{x * 1000:.1f}' if x < 1 else f'{x * 1000:.0f}'      # microseconds


def cls_label(c):
    return c.replace('+-', '±')


def params(cfg, rule):
    r = cfg.get('rel', {})
    if rule == 'SCA':
        return f"q {r['q']:g} T0 {r['T0']:g}"
    if rule == 'TEC':
        return f"Jv {r['jv']:+g} q {r['q']:g} T0 {r['T0']:g}"
    if rule == 'Onsager-kT':
        return f"κ {r['kappa']:g} q {r['q']:g} T0 {r['T0']:g}"
    if rule == 'Onsager-online':
        return f"λ {r['lam_k2000']:g} q {r['q']:g} T0 {r['T0']:g}"
    return ''


def per_instance(obj='E1'):
    key = 'tts_ms'
    L = [f"| Instance | Class | BKV | Target | " + ' | '.join(f'{SH[m]}: S*, p, TTS (µs)' for m in RULES) +
         ' | plain/kT | plain/online | TEC/kT | TEC/online |', '|' + '---|' * (4 + len(RULES) + 4)]
    for g, rec in INST.items():
        if 'ratios' not in rec:
            continue
        cells = []
        for m in RULES:
            c = rec['rules'][m][f'conf_{obj}']
            cells.append(f"{f(c['S'])}, {f(c.get('p_target'), 3)}, {ms(c[key])}")
        rs = []
        for k in ('plain/kT', 'plain/online', 'TEC/kT', 'TEC/online'):
            r = rec['ratios'][f'{obj}:{k}']
            ci = r.get('ci95')
            rs.append(f"{f(r['ratio'], 3)}" + (f" [{f(ci[0], 3)}, {f(ci[1], 3)}]" if ci else ''))
        L.append(f"| {g} | {cls_label(rec['cls'])} | {rec['BKV']} | {rec['target']} | " + ' | '.join(cells) + ' | ' + ' | '.join(rs) + ' |')
    return '\n'.join(L)


def configs_table():
    L = ['| Instance | ' + ' | '.join(f'{SH[m]} (E1 confirmed)' for m in RULES) + ' |', '|' + '---|' * (1 + len(RULES))]
    for g, rec in INST.items():
        if 'ratios' not in rec:
            continue
        cells = []
        for m in RULES:
            c = rec['rules'][m]['conf_E1']
            cells.append('–' if c['cfg'] is None else f"S {c['S']}: {params(c['cfg'], m)}")
        L.append(f'| {g} | ' + ' | '.join(cells) + ' |')
    return '\n'.join(L)


def summary_table():
    S = D['summary']
    L = ['| Ratio (TTS_ref / TTS_rule) | Estimator | Geometric mean | Min | Median | Max | n (both finite) | Rule faster (of n) | Only rule solves | Only ref solves | Neither |',
         '|---|---|---|---|---|---|---|---|---|---|---|']
    for obj in ('E1', 'E12'):
        for k in ('plain/online', 'plain/kT', 'plain/bestOnsager', 'TEC/online', 'TEC/kT', 'TEC/bestOnsager', 'plain/TEC'):
            v = S[f'{obj}:{k}']; s = v['stats']
            if s is None:
                L.append(f"| {k} | {obj} | – | – | – | – | 0 | – | {v['n_x_only']} | {v['n_ref_only']} | {v['n_neither']} |")
                continue
            L.append(f"| {k} | {'1 engine' if obj == 'E1' else '12-engine rounds'} | **{s['geomean']:.3f}** | {s['min']:.3f} | {s['median']:.3f} | "
                     f"{s['max']:.3f} | {s['n']} | {round(s['frac_gt1'] * s['n'])} | {v['n_x_only']} | {v['n_ref_only']} | {v['n_neither']} |")
    return '\n'.join(L)


def class_table(obj='E1'):
    S = D['summary']
    keys = ('plain/online', 'plain/kT', 'TEC/online', 'TEC/kT')
    L = ['| Class | n | ' + ' | '.join(f'{k}: wins, geomean' for k in keys) + ' |', '|' + '---|' * (2 + len(keys))]
    for c in CLASS_ORDER:
        row = []; n = None
        for k in keys:
            pc = S[f'{obj}:{k}']['per_class'].get(c)
            if pc is None:
                row.append('–'); continue
            n = pc['n']
            row.append(f"{pc['x_wins']}/{pc['n']}, {f(pc['stats']['geomean'], 3) if pc['stats'] else '–'}")
        L.append(f'| {cls_label(c)} | {n} | ' + ' | '.join(row) + ' |')
    return '\n'.join(L)


def mcs_table(methods, only_complete=True):
    L = ['| Instance | ' + ' | '.join(methods) + ' |', '|' + '---|' * (1 + len(methods))]
    for g, rec in INST.items():
        if only_complete and not all(m in rec['rules'] for m in methods):
            continue
        L.append(f'| {g} | ' + ' | '.join(f(rec['rules'][m]['mcs99_min'], 4) + (f" (S{rec['rules'][m]['mcs99_S']})" if rec['rules'][m]['mcs99_S'] else '')
                                        for m in methods) + ' |')
    return '\n'.join(L)


def cut_table(methods, S):
    L = [f'| Instance | BKV | ' + ' | '.join(methods) + ' |', '|' + '---|' * (2 + len(methods))]
    for g, rec in INST.items():
        if not all(m in rec['rules'] for m in methods):
            continue
        cells = []
        for m in methods:
            ps = rec['rules'][m]['per_S'][str(S)] if str(S) in rec['rules'][m]['per_S'] else rec['rules'][m]['per_S'][S]
            cells.append(f"{ps.get('mean_cut', float('nan')):.1f} / {f(ps.get('p_target'), 2)}")
        L.append(f"| {g} | {rec['BKV']} | " + ' | '.join(cells) + ' |')
    return '\n'.join(L)


def bkv_table(methods):
    L = ['| Instance | BKV | ' + ' | '.join(methods) + ' |', '|' + '---|' * (2 + len(methods))]
    for g, rec in INST.items():
        if not all(m in rec['rules'] for m in methods):
            continue
        L.append(f"| {g} | {rec['BKV']} | " + ' | '.join(f"{rec['rules'][m]['k_bkv_finals']} ({rec['rules'][m]['best_cut_finals']})"
                                                       for m in methods) + ' |')
    return '\n'.join(L)


GROUPS = (('Random, +1', ('R800+', 'R2000+', 'R1000+')), ('Random, ±1', ('R800+-', 'R2000+-')), ('Toroidal, ±1', ('T800+-', 'T2000+-')),
          ('Planar-like, +1', ('P800+', 'P2000+', 'P1000+')), ('Planar-like, ±1', ('P800+-', 'P2000+-')))


def group_stats(keys, obj, classes):
    vals = [rec['ratios'][f'{obj}:{keys}']['ratio'] for rec in INST.values() if 'ratios' in rec and rec['cls'] in classes]
    fin = [v for v in vals if v is not None and math.isfinite(v) and v > 0]
    wins = sum(1 for v in vals if v is not None and v > 1)
    gm = math.exp(sum(math.log(v) for v in fin) / len(fin)) if fin else None
    return len(vals), gm, wins, (min(fin) if fin else None), (max(fin) if fin else None)


def paper_table(obj='E12'):
    L = ['| G-set group (instances) | Onsager-kT vs plain SCA | Onsager-online vs plain SCA | Onsager-kT vs TEC | Onsager-online vs TEC |',
         '|---|---|---|---|---|']
    allc = tuple(c for _, cs in GROUPS for c in cs)
    for name, cs in GROUPS + (('All', allc),):
        cells = []; n = 0
        for k in ('plain/kT', 'plain/online', 'TEC/kT', 'TEC/online'):
            n, gm, w, lo, hi = group_stats(k, obj, cs)
            cells.append(f"{f(gm, 3)}× ({w}/{n} faster)")
        L.append(f'| {name} ({n}) | ' + ' | '.join(cells) + ' |')
    return '\n'.join(L)


if __name__ == '__main__':
    print('## paper table E12\n' + paper_table('E12') + '\n\n## paper table E1\n' + paper_table('E1'))
    print('## per-instance E1\n' + per_instance('E1') + '\n\n## per-instance E12\n' + per_instance('E12'))
    print('\n## configs\n' + configs_table())
    print('\n## summary\n' + summary_table())
    print('\n## classes E1\n' + class_table('E1') + '\n\n## classes E12\n' + class_table('E12'))
    print('\n## mcs99 rules\n' + mcs_table(RULES))
    print('\n## mcs99 all\n' + mcs_table(METHODS))
    print('\n## bkv all\n' + bkv_table(METHODS))
