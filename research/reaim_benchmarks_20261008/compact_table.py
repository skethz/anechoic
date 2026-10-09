"""Compact paper table (8 rows x 10 methods) from results_summary.json: per problem, mean normalized quality at the
ReAIM-matched budget (MCP S = 4000, GPP 4096, TSP 8192), feasibility (GPP, TSP), and steps to the 1% target
(geometric mean of the per-instance minimum MCS99 over solved instances, with the number solved).
Writes compact_table.md."""
import json
import math
from pathlib import Path

HERE = Path(__file__).resolve().parent
M = ('SA', 'SCA', 'TEC', 'APC-SCA', 'ReAIM ASA', 'aSB', 'bSB', 'dSB', 'Onsager-kT', 'Onsager-online')
H = ('SA', 'SCA', 'TEC', 'APC', 'ReAIM', 'aSB', 'bSB', 'dSB', 'Ons-κT', 'Ons-onl')
S = {'mcp': '4000', 'gpp': '4096', 'tsp': '8192'}


def geo(x):
    x = [v for v in x if math.isfinite(v)]
    return math.exp(sum(map(math.log, x)) / len(x)) if x else math.inf


def main():
    d = json.loads((HERE / 'results_summary.json').read_text())['summary']
    rows = []
    for prob, label in (('mcp', 'Max-Cut G1–G20 (20)'), ('gpp', 'GPP G1–5, G14–17 (9)'), ('tsp', 'TSP 6 TSPLIB (6)')):
        P = d[prob]; insts = list(P)
        q = [sum(P[i][m]['by_S'][S[prob]]['mean_quality'] for i in insts) / len(insts) for m in M]
        rows.append((f'{label}: quality', [f'{v:.3f}' for v in q]))
        if prob != 'mcp':
            f = [sum(P[i][m]['by_S'][S[prob]]['p_feasible'] for i in insts) / len(insts) for m in M]
            rows.append((f'{prob.upper()}: feasible', [f'{v:.2f}' for v in f]))
        st = []
        for m in M:
            vals = [float(P[i][m]['mcs99_min']) for i in insts]
            n = sum(1 for v in vals if math.isfinite(v)); g = geo(vals)
            st.append('–' if n == 0 else f'{g:,.0f} ({n}/{len(insts)})')
        rows.append((f'{prob.upper()}: steps to 1% target', st))
    L = ['| | ' + ' | '.join(H) + ' |', '|---|' + '---|' * len(H)]
    for name, vals in rows:
        L.append(f'| {name} | ' + ' | '.join(vals) + ' |')
    txt = '\n'.join(L)
    (HERE / 'compact_table.md').write_text(txt + '\n')
    print(txt)


if __name__ == '__main__':
    main()
