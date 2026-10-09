"""Addendum 4 / 4.1 Table 8b assembly: per method and problem, the source of every number and its label.
Sources: A3 finals (SA, STATICA sigma-transfer, bSB/dSB on Max-Cut and GPP, aSB dt 0.9, Onsager forms), A3.1 (APC-SCA,
ReAIM with their output rules), A4 (TEC sequential sigma-transfer, aSB dt 0.5 and stability rule), A4.1 (SB on TSP with
the ancillary spin). Writes table8b_a4.json and table8b_a4.md. Usage: python3 analyze_a4.py"""
import json
import math
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.dont_write_bytecode = True
sys.path.insert(0, str(HERE))
import grids_a2 as GA  # noqa: E402

S_REAIM = {'mcp': 4000, 'gpp': 4096, 'tsp': 8192}
NAME = {'mcp': 'Max-Cut G1–G20', 'gpp': 'GPP (9)', 'tsp': 'TSP (6)'}


def tag(prob, inst):
    return inst if prob == 'tsp' else f'G{inst}'


def src(entry, prob):
    """entry: (directory relative to HERE, method file name) or a dict per problem."""
    return entry[prob] if isinstance(entry, dict) else entry


ROWS = [  # (row label, method file, source per problem, label per problem)
    ('SA (Neal)', 'SA', 'results_a3/final', None),
    ('SCA (STATICA), σ-transfer [our rule]', 'SCA', 'results_a3/final', None),
    ('TEC sequential, σ-transfer [our rule]', 'TEC', 'results_a4/final/tec_seq', None),
    ('TEC synchronous [our misreading]', 'TEC', 'results_a3/final', None),
    ('APC-SCA (GPP: TSP column [our choice])', 'APC-SCA', 'results_a3/final_bv', None),
    ('ReAIM ASA (k set, ITER [our choice])', 'ReAIM_ASA', 'results_a3/final_bv', None),
    ('aSB Δt 0.9 M 2 (paper K2000 values)', 'aSB', {'mcp': 'results_a3/final', 'gpp': 'results_a3/final',
                                                    'tsp': 'results_a4/final/sb_anc/aSB_dt09'}, None),
    ('aSB Δt 0.5 M 2 [our choice]', 'aSB', {'mcp': 'results_a4/final/asb_dt05', 'gpp': 'results_a4/final/asb_dt05',
                                             'tsp': 'results_a4/final/sb_anc/aSB_dt05'}, None),
    ('bSB', 'bSB', {'mcp': 'results_a3/final', 'gpp': 'results_a3/final', 'tsp': 'results_a4/final/sb_anc/bSB'}, None),
    ('dSB', 'dSB', {'mcp': 'results_a3/final', 'gpp': 'results_a3/final', 'tsp': 'results_a4/final/sb_anc/dSB'}, None),
    ('Onsager-κT (ours)', 'Onsager-kT', 'results_a3/final', None),
    ('Onsager-online (ours)', 'Onsager-online', 'results_a3/final', None),
]
EXTRA = [('bSB, frozen ancilla [earlier choice]', 'bSB', 'results_a3/final', 'tsp'),
         ('dSB, frozen ancilla [earlier choice]', 'dSB', 'results_a3/final', 'tsp'),
         ('aSB stability rule [our rule]', 'aSB', 'results_a4/final/asb_stab', None)]


def load(dirpath, prob, inst, mfile, S):
    f = HERE / dirpath / prob / tag(prob, inst) / f'{mfile}_S{S}.json'
    return json.loads(f.read_text())['final'] if f.exists() else None


def stats_for(dirpath, mfile, prob):
    insts = GA.INSTANCES[prob]; S = S_REAIM[prob]
    q = []; f = []; mcs = {}; inv = 0; missing = 0
    for inst in insts:
        e = load(dirpath, prob, inst, mfile, S)
        if e is None:
            missing += 1; continue
        q.append(e['mean_quality']); f.append(e['p_feasible']); inv += e.get('n_invalid', 0)
        vals = [load(dirpath, prob, inst, mfile, s) for s in GA.S_LIST[prob]]
        mcs[tag(prob, inst)] = min(float(v['mcs99']) for v in vals if v is not None)
    if missing:
        return None
    solved = {t: v for t, v in mcs.items() if math.isfinite(v)}
    g = math.exp(np.mean([math.log(v) for v in solved.values()])) if solved else None
    return dict(quality=float(np.mean(q)), feasible=float(np.mean(f)), invalid_runs=inv, mcs=mcs, solved=len(solved),
                n=len(insts), geo_own=g)


def main():
    T = {}
    for label, mfile, s, _ in ROWS + [(a, b, c, None) for a, b, c, _ in EXTRA]:
        T[label] = {}
        for prob in ('mcp', 'gpp', 'tsp'):
            T[label][prob] = stats_for(src(s, prob), mfile, prob)
    # common solved Max-Cut set over the main rows that solve anything
    main_rows = [r[0] for r in ROWS]
    mcp_rows = [r for r in main_rows if T[r]['mcp'] and T[r]['mcp']['solved'] > 0]
    insts = [tag('mcp', g) for g in GA.MCP_INST]
    common = [t for t in insts if all(math.isfinite(T[r]['mcp']['mcs'][t]) for r in mcp_rows)]
    for r in main_rows:
        e = T[r]['mcp']
        e['geo_common'] = (math.exp(np.mean([math.log(e['mcs'][t]) for t in common]))
                           if common and r in mcp_rows else None)
    out = dict(table=T, mcp_common_set=common, mcp_rows_in_common=mcp_rows)
    (HERE / 'table8b_a4.json').write_text(json.dumps(out, indent=1, default=float) + '\n')
    f = lambda x, nd=3: '–' if x is None else f'{x:.{nd}f}'  # noqa: E731
    fs = lambda x: '–' if x is None or not math.isfinite(x) else f'{x:,.0f}'  # noqa: E731
    L = [f'Max-Cut common solved set ({len(common)} instances, over the rows that solve at least one): {common}\n',
         '| Method | Max-Cut quality | Max-Cut steps, common set | Max-Cut solved/20 | GPP quality | GPP feasible | GPP solved/9 | TSP quality | TSP valid | TSP solved/6 |',
         '|---|---|---|---|---|---|---|---|---|---|']
    for r in main_rows + [e[0] for e in EXTRA]:
        x = T[r]
        cells = []
        for prob in ('mcp', 'gpp', 'tsp'):
            e = x[prob]
            if e is None:
                cells += ['–'] * (3 if prob != 'mcp' else 3); continue
            if prob == 'mcp':
                cells += [f(e['quality']), fs(e.get('geo_common')), f"{e['solved']}/20"]
            else:
                cells += [f(e['quality']), f(e['feasible'], 2), f"{e['solved']}/{e['n']}"]
        L.append(f'| {r} | ' + ' | '.join(cells) + ' |')
    (HERE / 'table8b_a4.md').write_text('\n'.join(L) + '\n')
    print('\n'.join(L))


if __name__ == '__main__':
    main()
