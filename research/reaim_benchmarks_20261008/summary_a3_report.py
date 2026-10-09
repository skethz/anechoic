"""Report numbers for Addendum 3/3.1 (reads results_summary_a3b.json, results_summary_a3.json, results_summary.json and the
raw finals): the compact Table 8b, Max-Cut step comparisons on several common sets, A1 vs A3.1, and APC/ReAIM final-state
vs paper-output. Writes summary_a3_report.json and prints markdown."""
import json
import math
import sys
from pathlib import Path

import numpy as np

sys.dont_write_bytecode = True
HERE = Path(__file__).resolve().parent
M = ('SA', 'SCA', 'TEC', 'APC-SCA', 'ReAIM ASA', 'aSB', 'bSB', 'dSB', 'Onsager-kT', 'Onsager-online')
SH = ('SA (Neal)', 'SCA (STATICA)', 'TEC', 'APC-SCA', 'ReAIM ASA', 'aSB', 'bSB', 'dSB', 'Onsager-κT', 'Onsager-online')


def geo(xs):
    xs = [x for x in xs if math.isfinite(x) and x > 0]
    return math.exp(sum(map(math.log, xs)) / len(xs)) if xs else None


def fs(x):
    return '–' if x is None or not math.isfinite(x) else f'{x:,.0f}'


def main():
    B = json.loads((HERE / 'results_summary_a3b.json').read_text())['summary']
    F = json.loads((HERE / 'results_summary_a3.json').read_text())['summary']
    out = {}
    X = B['mcp']['rows']
    insts = list(X['SA']['mcs3'])
    solved = {m: [t for t in insts if math.isfinite(X[m]['mcs3'][t])] for m in M}
    sets = {
        'all ten methods': [t for t in insts if all(math.isfinite(X[m]['mcs3'][t]) for m in M)],
        'all methods except TEC (0/20) and aSB': [t for t in insts if all(math.isfinite(X[m]['mcs3'][t]) for m in M if m not in ('TEC', 'aSB'))],
        'SA, dSB, Onsager-online (C3)': [t for t in insts if all(math.isfinite(X[m]['mcs3'][t]) for m in ('SA', 'dSB', 'Onsager-online'))],
        'A1 common set (11)': ['G1', 'G2', 'G3', 'G4', 'G5', 'G6', 'G7', 'G8', 'G9', 'G10', 'G16'],
    }
    L = ['## Max-Cut steps to the 1% target (MCS99, min over S), geometric means over common solved sets\n']
    L.append('| Set (n) | ' + ' | '.join(SH) + ' |'); L.append('|---|' + '---|' * len(M))
    out['mcp_sets'] = {}
    for name, S in sets.items():
        row = []
        for m in M:
            v = [X[m]['mcs3'][t] for t in S]
            row.append(fs(geo(v)) if S and all(math.isfinite(x) for x in v) else ('–' if S else 'empty'))
        out['mcp_sets'][name] = dict(instances=S, geo={m: (geo([X[m]['mcs3'][t] for t in S]) if S and all(math.isfinite(X[m]['mcs3'][t]) for t in S) else None) for m in M})
        L.append(f'| {name} ({len(S)}) | ' + ' | '.join(row) + ' |')
    L.append('| solved / 20 | ' + ' | '.join(str(len(solved[m])) for m in M) + ' |')
    g = out['mcp_sets']
    c3 = g['SA, dSB, Onsager-online (C3)']['geo']; c11 = g['A1 common set (11)']['geo']
    out['C3'] = dict(online_over_SA_all20=c3['Onsager-online'] / c3['SA'], online_over_dSB_all20=c3['Onsager-online'] / c3['dSB'],
                     online_over_SA_A1set=(c11['Onsager-online'] / c11['SA']) if c11['SA'] else None,
                     online_over_dSB_A1set=(c11['Onsager-online'] / c11['dSB']) if c11['dSB'] else None,
                     kT_over_SA_all20=c3 and (geo([X['Onsager-kT']['mcs3'][t] for t in insts]) / c3['SA']),
                     kT_over_dSB_all20=geo([X['Onsager-kT']['mcs3'][t] for t in insts]) / c3['dSB'])
    L.append('\n## Max-Cut paired step ratios, method / Onsager form, over instances both solve (n; Onsager fewer steps on k)\n')
    L.append('| Method | vs Onsager-κT: geo ratio (n, k) | vs Onsager-online: geo ratio (n, k) |'); L.append('|---|---|---|')
    out['paired'] = {}
    for m, sh in zip(M, SH):
        cells = []
        for o in ('Onsager-kT', 'Onsager-online'):
            r = [X[m]['mcs3'][t] / X[o]['mcs3'][t] for t in insts if math.isfinite(X[m]['mcs3'][t]) and math.isfinite(X[o]['mcs3'][t])]
            gm = geo(r); k = sum(1 for x in r if x > 1)
            out['paired'][f'{m}/{o}'] = dict(geo=gm, n=len(r), onsager_fewer=k)
            cells.append(f"{gm:.2f}× ({len(r)}, {k})" if gm else f'– (0 solved by {sh})')
        L.append(f'| {sh} | ' + ' | '.join(cells) + ' |')
    L.append('\nC3: ' + json.dumps({k: round(v, 3) if isinstance(v, float) else v for k, v in out['C3'].items()}))
    # compact table
    K = ['| | ' + ' | '.join(SH) + ' |', '|---|' + '---|' * len(M)]
    for prob, S, lab in (('mcp', 4000, 'Max-Cut G1–G20: quality (S=4000)'), ('gpp', 4096, 'GPP (9): quality (S=4096)'),
                         ('tsp', 8192, 'TSP (6): quality (S=8192)')):
        R = B[prob]['rows']
        K.append(f'| {lab} | ' + ' | '.join(f"{R[m]['q3']:.3f}" for m in M) + ' |')
        if prob != 'mcp':
            K.append(f'| {prob.upper()}: feasible | ' + ' | '.join(f"{R[m]['f3']:.2f}" for m in M) + ' |')
        n = len(R['SA']['mcs3'])
        if prob == 'mcp':
            S3 = sets['SA, dSB, Onsager-online (C3)']
            K.append(f'| Max-Cut: steps, geo-mean over own solved (solved/20) | ' +
                     ' | '.join(f"{fs(R[m]['geo_own3'])} ({R[m]['solved3']}/20)" for m in M) + ' |')
        else:
            K.append(f'| {prob.upper()}: steps, own solved (solved/{n}) | ' +
                     ' | '.join(f"{fs(R[m]['geo_own3'])} ({R[m]['solved3']}/{n})" for m in M) + ' |')
    L.append('\n## Compact Table 8b (A3.1)\n'); L += K
    # A1 vs A3.1
    L.append('\n## A1 (tuned baselines) vs A3.1 (papers\' settings): quality (feasibility)\n')
    L.append('| Method | MCP A1 | MCP A3.1 | GPP A1 | GPP A3.1 | TSP A1 | TSP A3.1 |'); L.append('|---|---|---|---|---|---|---|')
    for m, sh in zip(M, SH):
        c = []
        for prob in ('mcp', 'gpp', 'tsp'):
            R = B[prob]['rows'][m]
            c.append(f"{R['q1']:.3f} ({R['f1']:.2f})"); c.append(f"{R['q3']:.3f} ({R['f3']:.2f})")
        L.append(f'| {sh} | ' + ' | '.join(c) + ' |')
    # output rule effect
    L.append('\n## APC-SCA and ReAIM: final state (A3, seed entry 41) vs the papers\' output rule (A3.1, seed entry 42)\n')
    L.append('| Method | Problem | final-state quality (feasible) | paper-output quality (feasible) |'); L.append('|---|---|---|---|')
    for m in ('APC-SCA', 'ReAIM ASA'):
        for prob in ('mcp', 'gpp', 'tsp'):
            a = F[prob]['rows'][m]; b = B[prob]['rows'][m]
            L.append(f"| {m} | {prob.upper()} | {a['q3']:.3f} ({a['f3']:.2f}) | {b['q3']:.3f} ({b['f3']:.2f}) |")
    txt = '\n'.join(L)
    (HERE / 'summary_a3_report.json').write_text(json.dumps(out, indent=1, default=str) + '\n')
    (HERE / 'summary_a3_report.md').write_text(txt + '\n')
    print(txt)


if __name__ == '__main__':
    main()
