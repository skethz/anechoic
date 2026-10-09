"""Addendum 2 analysis: Table 8b old (A1) vs new (A2), edge counts after re-tuning, extension logs, conclusions C1-C3.
Writes analysis_a2.txt and results_summary_a2.json. Usage: python3 analyze_a2.py"""
import json
import math
import sys
from pathlib import Path

import numpy as np

sys.dont_write_bytecode = True
HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import grids_a2 as GA  # noqa: E402

M = GA.METHODS
SHORT = {'SA': 'SA', 'SCA': 'SCA', 'TEC': 'TEC', 'APC-SCA': 'APC', 'ReAIM ASA': 'ReAIM', 'aSB': 'aSB', 'bSB': 'bSB',
         'dSB': 'dSB', 'Onsager-kT': 'Ons-kT', 'Onsager-online': 'Ons-onl'}
ENGINE = ('SCA', 'TEC', 'Onsager-kT', 'Onsager-online')
S_REAIM = {'mcp': 4000, 'gpp': 4096, 'tsp': 8192}
OUT = HERE / 'results_a2'


def tag(prob, inst):
    return inst if prob == 'tsp' else f'G{inst}'


def load_new():
    """{prob: {tag: {method: {S: final dict}}}} from results_a2/final."""
    out = {}
    for prob, insts in GA.INSTANCES.items():
        out[prob] = {}
        for inst in insts:
            t = tag(prob, inst); out[prob][t] = {}
            for m in M:
                out[prob][t][m] = {}
                for S in GA.S_LIST[prob]:
                    f = OUT / 'final' / prob / t / f"{m.replace(' ', '_')}_S{S}.json"
                    out[prob][t][m][S] = json.loads(f.read_text()) if f.exists() else None
    return out


def load_old():
    d = json.loads((HERE / 'results_summary.json').read_text())['summary']
    out = {}
    for prob in d:
        out[prob] = {}
        for t, ms in d[prob].items():
            out[prob][t] = {m: {int(S): e for S, e in ms[m]['by_S'].items()} for m in ms}
    return out


def geo(xs):
    xs = [x for x in xs if x is not None and math.isfinite(x) and x > 0]
    return math.exp(sum(math.log(x) for x in xs) / len(xs)) if xs else None


def table_rows(prob, old, new):
    insts = list(new[prob])
    S = S_REAIM[prob]
    rows = {}
    for m in M:
        q_o = [old[prob][t][m][S]['mean_quality'] for t in insts]
        f_o = [old[prob][t][m][S]['p_feasible'] for t in insts]
        q_n = [new[prob][t][m][S]['final']['mean_quality'] for t in insts]
        f_n = [new[prob][t][m][S]['final']['p_feasible'] for t in insts]
        mcs_o = [min(float(old[prob][t][m][s]['mcs99']) for s in GA.S_LIST[prob]) for t in insts]
        mcs_n = [min(float(new[prob][t][m][s]['final']['mcs99']) for s in GA.S_LIST[prob]) for t in insts]
        rows[m] = dict(q_old=float(np.mean(q_o)), q_new=float(np.mean(q_n)), feas_old=float(np.mean(f_o)),
                       feas_new=float(np.mean(f_n)), feas_new_range=(float(min(f_n)), float(max(f_n))),
                       feas_old_range=(float(min(f_o)), float(max(f_o))),
                       mcs_old=dict(zip(insts, mcs_o)), mcs_new=dict(zip(insts, mcs_n)),
                       solved_old=sum(1 for x in mcs_o if math.isfinite(x)), solved_new=sum(1 for x in mcs_n if math.isfinite(x)),
                       q_new_per_inst=dict(zip(insts, q_n)), feas_new_per_inst=dict(zip(insts, f_n)),
                       q_old_per_inst=dict(zip(insts, q_o)))
    common_o = [t for t in insts if all(math.isfinite(rows[m]['mcs_old'][t]) for m in M)]
    common_n = [t for t in insts if all(math.isfinite(rows[m]['mcs_new'][t]) for m in M)]
    for m in M:
        rows[m]['geo_common_old'] = geo([rows[m]['mcs_old'][t] for t in common_o]) if common_o else None
        rows[m]['geo_common_new'] = geo([rows[m]['mcs_new'][t] for t in common_n]) if common_n else None
        rows[m]['geo_own_old'] = geo(list(rows[m]['mcs_old'].values()))
        rows[m]['geo_own_new'] = geo(list(rows[m]['mcs_new'].values()))
    return rows, common_o, common_n


def fmt(x, nd=3):
    return '–' if x is None else f'{x:.{nd}f}'


def fmt_s(x):
    return '–' if x is None or not math.isfinite(x) else f'{x:,.0f}'


def edges_after(R=GA.MAX_ROUNDS):
    """Edge counts of the final A2 selections, and the extension logs (grids_r<R>.json written by run_a2)."""
    import run_a2 as RA
    grids, report = RA.build_grids(R)
    res = {}
    for (fam, m), g in grids.items():
        sels = []
        for prob, inst in GA.family_instances(fam):
            for S in GA.S_LIST[prob]:
                idx, rec, tied = RA.a2_selection(prob, inst, m, S, R, g)
                sels.append((g.points[g.order[idx]], tied))
        c = GA.edge_counts(g, sels)
        res[f'{fam}|{m}'] = dict(size=g.size(), cap=g.cap(), a1=g.a1_size, counts=c, log=g.log,
                                 residual=[(k, e, v[f'frac_{e}']) for k, v in c.items() for e in ('min', 'max') if v[f'frac_{e}'] > GA.THRESH],
                                 n_tied=sum(1 for _, t in sels if t))
    return res


def selection_changes(new):
    """Engine-rule precision-study schedules (S_exp rule of precision.py) under A1 and A2."""
    import precision as PR
    out = {}
    old_s = {(s['problem'], str(s['instance']), s['rule']): s for s in json.loads((HERE / 'precision' / 'schedules.json').read_text())}
    for prob in ('mcp', 'gpp', 'tsp'):
        for t in new[prob]:
            inst = t if prob == 'tsp' else int(t[1:])
            for r in PR.RULES:
                recs = [(S, new[prob][t][r][S]['final'] | dict(mean_quality=new[prob][t][r][S]['final']['mean_quality']),
                         new[prob][t][r][S]['selected_cfg']) for S in GA.S_LIST[prob]]
                S, c, why = PR.pick(recs)
                o = old_s[(prob, str(inst), r)]
                changed = (S != o['S']) or any(c.get(k) != o['cfg'].get(k) for k in ('q', 'T0', 'tfin', 'jv', 'kappa', 'lam'))
                out[f'{prob}|{t}|{r}'] = dict(S_old=o['S'], S_new=S, why=why, changed=bool(changed), cfg_new=c, cfg_old=o['cfg'])
    return out


def main():
    old, new = load_old(), load_new()
    L = ['# Addendum 2: re-tuned results (A2) against the first run (A1)\n']
    summ = {}
    for prob in ('mcp', 'gpp', 'tsp'):
        rows, co, cn = table_rows(prob, old, new)
        summ[prob] = dict(rows=rows, common_old=co, common_new=cn)
        S = S_REAIM[prob]
        L.append(f'\n## {prob.upper()}: quality and feasibility at S = {S}; steps (MCS99, min over S)\n')
        L.append(f'Common solved set: A1 {len(co)} instances {co}; A2 {len(cn)} instances {cn}\n')
        L.append('| Method | quality A1 | quality A2 | feasible A1 | feasible A2 | solved A1 | solved A2 | steps A1 (common) | '
                 'steps A2 (common) | steps A2 (own solved) |')
        L.append('|---|---|---|---|---|---|---|---|---|---|')
        n = len(new[prob])
        for m in M:
            r = rows[m]
            L.append(f"| {m} | {r['q_old']:.4f} | {r['q_new']:.4f} | {r['feas_old']:.3f} | {r['feas_new']:.3f} | "
                     f"{r['solved_old']}/{n} | {r['solved_new']}/{n} | {fmt_s(r['geo_common_old'])} | {fmt_s(r['geo_common_new'])} | "
                     f"{fmt_s(r['geo_own_new'])} |")
        L.append(f'\nPer instance, A2, S = {S}: mean quality (feasible fraction)\n')
        L.append('| Inst | ' + ' | '.join(SHORT[m] for m in M) + ' |')
        L.append('|---|' + '---|' * len(M))
        for t in new[prob]:
            L.append(f'| {t} | ' + ' | '.join(f"{rows[m]['q_new_per_inst'][t]:.3f} ({rows[m]['feas_new_per_inst'][t]:.2f})" for m in M) + ' |')
        L.append('\nPer instance, A2: MCS99 (min over S)\n')
        L.append('| Inst | ' + ' | '.join(SHORT[m] for m in M) + ' |')
        L.append('|---|' + '---|' * len(M))
        for t in new[prob]:
            L.append(f'| {t} | ' + ' | '.join(fmt_s(rows[m]['mcs_new'][t]) for m in M) + ' |')
        if prob == 'tsp':
            L.append('\nTSP ARPD over feasible runs, S = 8192, A2 (A1 in parentheses)\n')
            L.append('| Inst | ' + ' | '.join(SHORT[m] for m in M) + ' |')
            L.append('|---|' + '---|' * len(M))
            for t in new[prob]:
                cells = []
                for m in M:
                    a2 = new[prob][t][m][8192]['final'].get('arpd_feasible'); a1 = old[prob][t][m][8192].get('arpd_feasible')
                    cells.append(f'{fmt(a2, 1)} ({fmt(a1, 1)})')
                L.append(f'| {t} | ' + ' | '.join(cells) + ' |')
    # C1-C3
    G = summ['gpp']['rows']; T = summ['tsp']['rows']; X = summ['mcp']['rows']
    C = {}
    C['C1_engine_gpp_feasible_A2'] = {m: (G[m]['feas_new'], G[m]['feas_new_range']) for m in ENGINE}
    C['C1_engine_gpp_feasible_A1'] = {m: (G[m]['feas_old'], G[m]['feas_old_range']) for m in ENGINE}
    C['C2_engine_tsp_quality_A2'] = {m: T[m]['q_new'] for m in ENGINE}
    C['C2_tsp_best_method_per_instance_A2'] = {t: max(M, key=lambda m: T[m]['q_new_per_inst'][t]) for t in summ['tsp']['rows']['SA']['q_new_per_inst']}
    C['C2_tsp_quality_A2'] = {m: T[m]['q_new'] for m in M}
    on = X['Onsager-online']
    C['C3_common_set_A1'] = summ['mcp']['common_old']; C['C3_common_set_A2'] = summ['mcp']['common_new']
    C['C3_ratio_online_over_SA_A1'] = on['geo_common_old'] / X['SA']['geo_common_old']
    C['C3_ratio_online_over_dSB_A1'] = on['geo_common_old'] / X['dSB']['geo_common_old']
    if on['geo_common_new']:
        C['C3_ratio_online_over_SA_A2'] = on['geo_common_new'] / X['SA']['geo_common_new']
        C['C3_ratio_online_over_dSB_A2'] = on['geo_common_new'] / X['dSB']['geo_common_new']
        C['C3_ratio_online_over_kT_A2'] = on['geo_common_new'] / X['Onsager-kT']['geo_common_new']
        C['C3_ratio_plain_over_kT_A2'] = X['SCA']['geo_common_new'] / X['Onsager-kT']['geo_common_new']
        C['C3_ratio_plain_over_online_A2'] = X['SCA']['geo_common_new'] / on['geo_common_new']
    # paired per-instance ratios over all 20 (both finite)
    def paired(a, b, key):
        r = [X[a][key][t] / X[b][key][t] for t in X[a][key] if math.isfinite(X[a][key][t]) and math.isfinite(X[b][key][t])]
        return dict(geo=geo(r), n=len(r), b_fewer=sum(1 for x in r if x > 1))
    for a, b in (('SCA', 'Onsager-kT'), ('SCA', 'Onsager-online'), ('TEC', 'Onsager-kT'), ('Onsager-online', 'SA'),
                 ('Onsager-online', 'dSB'), ('Onsager-kT', 'SA'), ('Onsager-kT', 'dSB')):
        C[f'paired_{a}_over_{b}_A1'] = paired(a, b, 'mcs_old'); C[f'paired_{a}_over_{b}_A2'] = paired(a, b, 'mcs_new')
    L.append('\n## Conclusions C1-C3 (C4 from the precision re-run)\n')
    L.append('```\n' + json.dumps(C, indent=1, default=str) + '\n```')
    E = edges_after()
    L.append('\n## Edge counts after re-tuning (final A2 selections; > 20% flagged)\n')
    for k, v in E.items():
        fam, m = k.split('|')
        cs = '; '.join(f"{a}: {c['min']}/{c['n']} min, {c['max']}/{c['n']} max" for a, c in v['counts'].items())
        flag = (' RESIDUAL: ' + ', '.join(f'{a}-{e} {f:.0%}' for a, e, f in v['residual'])) if v['residual'] else ''
        L.append(f"- {fam} {m} (grid {v['size']}/{v['cap']}, A1 {v['a1']}, tied selections {v['n_tied']}): {cs}{flag}")
    L.append('\n## Extension logs\n')
    for k, v in E.items():
        if v['log']:
            L.append(f"- {k}: " + '; '.join(f"r{e['round']} {e['axis']}-{e['end']}={e.get('value')} ({e['result']})" for e in v['log']))
    SC = selection_changes(new)
    L.append('\n## Engine-rule precision-study schedules: A1 vs A2 (changed?)\n')
    for prob in ('mcp', 'gpp', 'tsp'):
        ch = [k for k, v in SC.items() if k.startswith(prob) and v['changed']]
        L.append(f'- {prob}: {len(ch)} of {sum(1 for k in SC if k.startswith(prob))} changed')
    txt = '\n'.join(L)
    (HERE / 'analysis_a2.txt').write_text(txt + '\n')
    (HERE / 'results_summary_a2.json').write_text(json.dumps(dict(summary=summ, conclusions=C, edges=E, schedules=SC),
                                                             indent=1, default=str) + '\n')
    print(txt)


if __name__ == '__main__':
    main()
