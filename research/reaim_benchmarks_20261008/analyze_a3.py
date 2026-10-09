"""Addendum 3 analysis: Table 8b with the baselines at their papers' settings (A3) against the first run (A1), the
conclusions C1-C4, and the comparison with ReAIM's own Tables IV-VI. Writes analysis_a3.txt, results_summary_a3.json,
compact_table_a3.md. Usage: python3 analyze_a3.py"""
import json
import math
import sys
from pathlib import Path

import numpy as np

sys.dont_write_bytecode = True
HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import grids_a2 as GA  # noqa: E402  (instance lists, budgets)

M = GA.METHODS
SHORT = {'SA': 'SA(Neal)', 'SCA': 'SCA(STATICA)', 'TEC': 'TEC', 'APC-SCA': 'APC', 'ReAIM ASA': 'ReAIM', 'aSB': 'aSB',
         'bSB': 'bSB', 'dSB': 'dSB', 'Onsager-kT': 'Ons-kT', 'Onsager-online': 'Ons-onl'}
ENGINE = ('SCA', 'TEC', 'Onsager-kT', 'Onsager-online')
S_REAIM = {'mcp': 4000, 'gpp': 4096, 'tsp': 8192}
BKV = json.loads((HERE.parent / 'gset_20261007' / 'bkv.json').read_text())['instances']
REF = json.loads((HERE / 'data' / 'gpp_reference.json').read_text())['instances']
TSPOPT = {'gr17': 2085, 'gr21': 2707, 'gr24': 1272, 'fri26': 937, 'bayg29': 1610, 'bays29': 2020}
COLS = ('SFG', 'MFG', 'SFA', 'MFA', 'ASF', 'AMF', 'ASA', 'Neal', 'Tabu')
REAIM_T4 = {  # ReAIM Table IV, best of 20 cut values; last entry = Best [7]
    'G1': (11552, 11613, 11561, 11620, 11550, 11624, 11624, 11624, 11589, 11624),
    'G2': (11566, 11592, 11561, 11617, 11571, 11616, 11615, 11620, 11600, 11620),
    'G6': (2145, 2167, 2123, 2176, 2137, 2178, 2168, 2177, 2134, 2178),
    'G7': (1967, 1991, 1939, 1996, 1985, 1999, 2006, 2006, 1984, 2006),
    'G10': (1941, 1990, 1954, 1998, 1970, 1990, 1991, 1996, 1958, 2000),
    'G11': (358, 452, 550, 558, 528, 530, 554, 560, 550, 564),
    'G12': (328, 456, 546, 550, 532, 530, 550, 554, 540, 556),
    'G13': (308, 458, 564, 576, 558, 548, 576, 576, 564, 582),
    'G14': (2926, 2994, 3022, 3045, 3015, 3040, 3046, 3056, 3003, 3064),
    'G19': (827, 884, 875, 895, 860, 899, 900, 893, 899, 906),
    'G20': (837, 925, 895, 937, 908, 936, 939, 940, 936, 941)}
REAIM_T5 = {  # ReAIM Table V, best of 20 cut values (lower is better)
    'G1': (7836, 8466, 7665, 7650, 7696, 7683, 7645, 7728, 7597), 'G2': (7799, 8503, 7676, 7672, 7659, 7665, 7611, 7763, 7616),
    'G3': (7830, 8467, 7656, 7649, 7657, 7693, 7588, 7718, 7626), 'G4': (7849, 8481, 7668, 7683, 7673, 7651, 7635, 7707, 7620),
    'G5': (7806, 8471, 7671, 7683, 7676, 7697, 7664, 7751, 7612), 'G14': (1318, 2065, 1144, 1136, 1143, 1142, 1115, 1251, 1150),
    'G15': (1373, 2059, 1131, 1149, 1148, 1121, 1126, 1284, 1188), 'G16': (1372, 2074, 1143, 1114, 1143, 1116, 1085, 1254, 1125),
    'G17': (1306, 2043, 1123, 1115, 1140, 1126, 1075, 1226, 1121)}
REAIM_T6 = {  # ReAIM Table VI, ARPD (%) over 20 runs
    'bayg29': (31.44, 92.73, 91.50, 43.48, 28.97, 39.93, 31.00, 26.27, 29.88),
    'bays29': (21.31, 83.66, 98.53, 49.52, 25.38, 38.46, 28.51, 31.29, 36.54),
    'fri26': (38.54, 103.00, 90.48, 48.94, 36.56, 40.37, 43.99, 29.35, 30.52),
    'gr17': (10.28, 58.98, 50.67, 21.60, 10.44, 16.41, 16.38, 14.92, 44.84),
    'gr21': (36.58, 108.83, 78.71, 39.59, 25.12, 32.40, 32.09, 32.95, 10.68),
    'gr24': (17.22, 71.25, 80.35, 40.79, 19.00, 26.14, 21.84, 30.82, 23.11)}


def tag(prob, inst):
    return inst if prob == 'tsp' else f'G{inst}'


def load_a3():
    out = {}
    for prob, insts in GA.INSTANCES.items():
        out[prob] = {}
        for inst in insts:
            t = tag(prob, inst); out[prob][t] = {}
            for m in M:
                out[prob][t][m] = {}
                for S in GA.S_LIST[prob]:
                    f = HERE / 'results_a3' / 'final' / prob / t / f"{m.replace(' ', '_')}_S{S}.json"
                    out[prob][t][m][S] = json.loads(f.read_text()) if f.exists() else None
    return out


def load_a1():
    d = json.loads((HERE / 'results_summary.json').read_text())['summary']
    return {prob: {t: {m: {int(S): e for S, e in ms[m]['by_S'].items()} for m in ms} for t, ms in d[prob].items()} for prob in d}


def geo(xs):
    xs = [x for x in xs if x is not None and math.isfinite(x) and x > 0]
    return math.exp(sum(math.log(x) for x in xs) / len(xs)) if xs else None


def fs(x):
    return '–' if x is None or not math.isfinite(x) else f'{x:,.0f}'


def rows(prob, a1, a3):
    insts = list(a3[prob]); S = S_REAIM[prob]
    R = {}
    for m in M:
        e3 = [a3[prob][t][m][S]['final'] for t in insts]
        e1 = [a1[prob][t][m][S] for t in insts]
        mcs3 = {t: min(float(a3[prob][t][m][s]['final']['mcs99']) for s in GA.S_LIST[prob]) for t in insts}
        mcs1 = {t: min(float(a1[prob][t][m][s]['mcs99']) for s in GA.S_LIST[prob]) for t in insts}
        R[m] = dict(q3=float(np.mean([e['mean_quality'] for e in e3])), f3=float(np.mean([e['p_feasible'] for e in e3])),
                    q1=float(np.mean([e['mean_quality'] for e in e1])), f1=float(np.mean([e['p_feasible'] for e in e1])),
                    f3_range=(float(min(e['p_feasible'] for e in e3)), float(max(e['p_feasible'] for e in e3))),
                    q3_inst={t: e['mean_quality'] for t, e in zip(insts, e3)}, f3_inst={t: e['p_feasible'] for t, e in zip(insts, e3)},
                    mcs3=mcs3, mcs1=mcs1, solved3=sum(1 for v in mcs3.values() if math.isfinite(v)),
                    solved1=sum(1 for v in mcs1.values() if math.isfinite(v)),
                    invalid3=sum(e.get('n_invalid', 0) for e in e3))
    common3 = [t for t in insts if all(math.isfinite(R[m]['mcs3'][t]) for m in M)]
    common1 = [t for t in insts if all(math.isfinite(R[m]['mcs1'][t]) for m in M)]
    for m in M:
        R[m]['geo_common3'] = geo([R[m]['mcs3'][t] for t in common3]) if common3 else None
        R[m]['geo_common1'] = geo([R[m]['mcs1'][t] for t in common1]) if common1 else None
        R[m]['geo_own3'] = geo(list(R[m]['mcs3'].values()))
    return R, common3, common1


def best_of_20(vals, feas, prob):
    """Median over the 12 disjoint blocks of 20 runs (runs 0-239) of the block's best feasible value; and blocks w/o feasible."""
    best = []; nofe = 0
    for k in range(12):
        v = np.array(vals[20 * k:20 * k + 20]); f = np.array(feas[20 * k:20 * k + 20], bool)
        if not f.any():
            nofe += 1; continue
        best.append(v[f].max() if prob == 'mcp' else v[f].min())
    return (float(np.median(best)) if best else None), nofe


def arpd_blocks(vals, feas, opt):
    """ARPD of each 20-run block over its feasible runs; median over blocks; total feasible fraction."""
    a = []
    for k in range(12):
        v = np.array(vals[20 * k:20 * k + 20], float); f = np.array(feas[20 * k:20 * k + 20], bool)
        if f.any():
            a.append(float((100.0 * (v[f] - opt) / opt).mean()))
    return (float(np.median(a)) if a else None)


def main():
    a1, a3 = load_a1(), load_a3()
    L = ['# Addendum 3: Table 8b with the baselines at their papers\' own settings (A3) vs the first run (A1)\n']
    summ = {}
    for prob in ('mcp', 'gpp', 'tsp'):
        R, c3, c1 = rows(prob, a1, a3)
        summ[prob] = dict(rows=R, common3=c3, common1=c1)
        n = len(a3[prob]); S = S_REAIM[prob]
        L.append(f'\n## {prob.upper()} (quality/feasibility at S = {S}; steps = MCS99 min over S)\n')
        L.append(f'Common solved set: A3 {len(c3)} {c3}; A1 {len(c1)} {c1}\n')
        L.append('| Method | quality A3 | quality A1 | feasible A3 | feasible A1 | solved A3 | solved A1 | steps A3 (common) | '
                 'steps A1 (common A1) | steps A3 (own solved) | invalid runs A3 |')
        L.append('|---|---|---|---|---|---|---|---|---|---|---|')
        for m in M:
            r = R[m]
            L.append(f"| {m} | {r['q3']:.4f} | {r['q1']:.4f} | {r['f3']:.3f} | {r['f1']:.3f} | {r['solved3']}/{n} | {r['solved1']}/{n} | "
                     f"{fs(r['geo_common3'])} | {fs(r['geo_common1'])} | {fs(r['geo_own3'])} | {r['invalid3']} |")
        L.append(f'\nPer instance (A3, S = {S}): mean quality (feasible)\n')
        L.append('| Inst | ' + ' | '.join(SHORT[m] for m in M) + ' |'); L.append('|---|' + '---|' * len(M))
        for t in a3[prob]:
            L.append(f'| {t} | ' + ' | '.join(f"{R[m]['q3_inst'][t]:.3f} ({R[m]['f3_inst'][t]:.2f})" for m in M) + ' |')
        L.append('\nPer instance (A3): MCS99 min over S\n')
        L.append('| Inst | ' + ' | '.join(SHORT[m] for m in M) + ' |'); L.append('|---|' + '---|' * len(M))
        for t in a3[prob]:
            L.append(f'| {t} | ' + ' | '.join(fs(R[m]['mcs3'][t]) for m in M) + ' |')
    # dt selections of bSB/dSB
    L.append('\n## bSB/dSB: dt selected by the paper\'s procedure (counts over instances x budgets)\n')
    for prob in ('mcp', 'gpp', 'tsp'):
        for m in ('bSB', 'dSB'):
            from collections import Counter
            c = Counter(a3[prob][t][m][S]['cfg']['dt'] for t in a3[prob] for S in GA.S_LIST[prob])
            L.append(f'- {prob} {m}: ' + ', '.join(f'dt={k}: {v}' for k, v in sorted(c.items())))
    # conclusions
    X, G, T = summ['mcp']['rows'], summ['gpp']['rows'], summ['tsp']['rows']
    C = {}
    C['C1_gpp_feasible_A3'] = {m: (G[m]['f3'], G[m]['f3_range']) for m in M}
    C['C2_tsp_quality_A3'] = {m: T[m]['q3'] for m in M}
    C['C2_tsp_best_per_instance_A3'] = {t: max(M, key=lambda m: T[m]['q3_inst'][t]) for t in T['SA']['q3_inst']}
    on = X['Onsager-online']['geo_common3']
    C['C3_common_set_A3'] = summ['mcp']['common3']
    if on:
        for b in M:
            C[f'C3_online_over_{b}_A3_common'] = on / X[b]['geo_common3']
    def paired(a, b):
        r = [X[a]['mcs3'][t] / X[b]['mcs3'][t] for t in X[a]['mcs3'] if math.isfinite(X[a]['mcs3'][t]) and math.isfinite(X[b]['mcs3'][t])]
        return dict(geo=geo(r), n=len(r), b_fewer_steps=sum(1 for x in r if x > 1))
    for a, b in (('Onsager-online', 'SA'), ('Onsager-online', 'dSB'), ('Onsager-kT', 'SA'), ('Onsager-kT', 'dSB'),
                 ('SCA', 'Onsager-kT'), ('SCA', 'Onsager-online'), ('TEC', 'Onsager-kT'), ('APC-SCA', 'Onsager-kT'),
                 ('ReAIM ASA', 'Onsager-kT'), ('bSB', 'Onsager-kT'), ('aSB', 'Onsager-kT')):
        C[f'paired_{a}_over_{b}'] = paired(a, b)
    L.append('\n## Conclusions (numbers)\n\n```\n' + json.dumps(C, indent=1, default=str) + '\n```')
    # ReAIM comparison
    L.append('\n## Comparison with ReAIM\'s own tables (our noise-free runs; ReAIM includes ReRAM noise)\n')
    L.append('Max-Cut (Table IV, best of 20 cut values, 4096 iterations; ours: median over 12 blocks of 20 runs at S = 4000)\n')
    L.append('| Inst | ours ReAIM ASA | ReAIM ASA (paper) | ours SA (Neal) | Neal (paper) | BKV |'); L.append('|---|---|---|---|---|---|')
    cmpd = {}
    for t, row in REAIM_T4.items():
        r3 = a3['mcp'][t]['ReAIM ASA'][4000]; n3 = a3['mcp'][t]['SA'][4000]
        br = best_of_20(r3['values'], r3['feasible_runs'], 'mcp')[0]; bn = best_of_20(n3['values'], n3['feasible_runs'], 'mcp')[0]
        cmpd[f'mcp/{t}'] = dict(ours_reaim=br, paper_asa=row[6], ours_neal=bn, paper_neal=row[7], bkv=row[9])
        L.append(f'| {t} | {br:.0f} | {row[6]} | {bn:.0f} | {row[7]} | {row[9]} |')
    L.append('\nGPP (Table V, best of 20 cut values, lower is better; ours: median over 12 blocks of the best balanced cut at S = 4096)\n')
    L.append('| Inst | ours ReAIM ASA (blocks w/o feasible) | ReAIM ASA (paper) | ours SA (Neal) | Neal (paper) | our reference R |')
    L.append('|---|---|---|---|---|---|')
    for t, row in REAIM_T5.items():
        r3 = a3['gpp'][t]['ReAIM ASA'][4096]; n3 = a3['gpp'][t]['SA'][4096]
        br, nr = best_of_20(r3['values'], r3['feasible_runs'], 'gpp'); bn, nn = best_of_20(n3['values'], n3['feasible_runs'], 'gpp')
        cmpd[f'gpp/{t}'] = dict(ours_reaim=br, ours_reaim_blocks_without_feasible=nr, paper_asa=row[6], ours_neal=bn,
                                ours_neal_blocks_without_feasible=nn, paper_neal=row[7], R=REF[t]['R'])
        L.append(f"| {t} | {('%.0f' % br) if br is not None else '–'} ({nr}) | {row[6]} | {('%.0f' % bn) if bn is not None else '–'} ({nn}) | {row[7]} | {REF[t]['R']} |")
    L.append('\nTSP (Table VI, ARPD %, 8192 iterations; ours: median over 12 blocks of 20 runs of the block ARPD over valid tours, '
             'and the valid-tour fraction)\n')
    L.append('| Inst | ours ReAIM ASA (valid) | ReAIM ASA (paper) | ours SA (Neal) (valid) | Neal (paper) |'); L.append('|---|---|---|---|---|')
    for t, row in REAIM_T6.items():
        r3 = a3['tsp'][t]['ReAIM ASA'][8192]; n3 = a3['tsp'][t]['SA'][8192]
        ar = arpd_blocks(r3['values'], r3['feasible_runs'], TSPOPT[t]); an = arpd_blocks(n3['values'], n3['feasible_runs'], TSPOPT[t])
        cmpd[f'tsp/{t}'] = dict(ours_reaim=ar, ours_reaim_valid=r3['final']['p_feasible'], paper_asa=row[6], ours_neal=an,
                                ours_neal_valid=n3['final']['p_feasible'], paper_neal=row[7])
        f = lambda x: '–' if x is None else f'{x:.1f}'  # noqa: E731
        L.append(f"| {t} | {f(ar)} ({r3['final']['p_feasible']:.2f}) | {row[6]} | {f(an)} ({n3['final']['p_feasible']:.2f}) | {row[7]} |")
    # compact table
    K = ['| | ' + ' | '.join(SHORT[m] for m in M) + ' |', '|---|' + '---|' * len(M)]
    K.append('| Max-Cut G1–G20: quality (S=4000) | ' + ' | '.join(f"{X[m]['q3']:.3f}" for m in M) + ' |')
    cs = summ['mcp']['common3']
    K.append(f'| Max-Cut: steps to 1% target, geo-mean over {len(cs)} common (solved/20) | ' +
             ' | '.join(f"{fs(X[m]['geo_common3'])} ({X[m]['solved3']}/20)" for m in M) + ' |')
    K.append('| GPP (9): quality (S=4096) | ' + ' | '.join(f"{G[m]['q3']:.3f}" for m in M) + ' |')
    K.append('| GPP: feasible | ' + ' | '.join(f"{G[m]['f3']:.2f}" for m in M) + ' |')
    K.append('| GPP: steps (own solved; solved/9) | ' + ' | '.join(f"{fs(G[m]['geo_own3'])} ({G[m]['solved3']}/9)" for m in M) + ' |')
    K.append('| TSP (6): quality (S=8192) | ' + ' | '.join(f"{T[m]['q3']:.3f}" for m in M) + ' |')
    K.append('| TSP: feasible | ' + ' | '.join(f"{T[m]['f3']:.2f}" for m in M) + ' |')
    K.append('| TSP: steps (own solved; solved/6) | ' + ' | '.join(f"{fs(T[m]['geo_own3'])} ({T[m]['solved3']}/6)" for m in M) + ' |')
    (HERE / 'compact_table_a3.md').write_text('\n'.join(K) + '\n')
    L.append('\n## Compact table (A3)\n\n' + '\n'.join(K))
    txt = '\n'.join(L)
    (HERE / 'analysis_a3.txt').write_text(txt + '\n')
    (HERE / 'results_summary_a3.json').write_text(json.dumps(dict(summary=summ, conclusions=C, reaim_comparison=cmpd),
                                                             indent=1, default=str) + '\n')
    print(txt)


if __name__ == '__main__':
    main()
