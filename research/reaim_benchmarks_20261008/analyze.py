"""Analysis of the held-out finals (PROTOCOL.md measures). Writes results_summary.json and analysis.txt (markdown tables).
MCP: four engine rules from research/gset_20261007 (pre-registered grids); SA, APC-SCA, ReAIM, aSB, bSB, dSB from its
Part B (G1, G6, G11, G14, G18) and from results/mcp_ext (the other 15 of G1-G20)."""
import json
import math
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
GSET = HERE.parent / 'gset_20261007'
METHODS = ('SA', 'SCA', 'TEC', 'APC-SCA', 'ReAIM ASA', 'aSB', 'bSB', 'dSB', 'Onsager-kT', 'Onsager-online')
SHORT = {'SA': 'SA', 'SCA': 'SCA', 'TEC': 'TEC', 'APC-SCA': 'APC', 'ReAIM ASA': 'ReAIM', 'aSB': 'aSB', 'bSB': 'bSB',
         'dSB': 'dSB', 'Onsager-kT': 'Ons-kT', 'Onsager-online': 'Ons-onl'}
ENGINE = ('SCA', 'TEC', 'Onsager-kT', 'Onsager-online')
S_LIST = {'mcp': (250, 500, 1000, 2000, 4000), 'gpp': (256, 512, 1024, 2048, 4096), 'tsp': (512, 1024, 2048, 4096, 8192)}
S_REAIM = {'mcp': 4000, 'gpp': 4096, 'tsp': 8192}
MCP_INST = tuple(range(1, 21)); SUBSET = (1, 6, 11, 14, 18)
GPP_INST = (1, 2, 3, 4, 5, 14, 15, 16, 17)
TSP_INST = ('gr17', 'gr21', 'gr24', 'fri26', 'bayg29', 'bays29')
BKV = json.loads((GSET / 'bkv.json').read_text())['instances']
REF = json.loads((HERE / 'data' / 'gpp_reference.json').read_text())['instances']
TSPOPT = {'gr17': 2085, 'gr21': 2707, 'gr24': 1272, 'fri26': 937, 'bayg29': 1610, 'bays29': 2020}
REAIM_T6 = {  # ReAIM Table VI ARPD (%), columns SFG MFG SFA MFA ASF AMF ASA Neal Tabu
    'bayg29': (31.44, 92.73, 91.50, 43.48, 28.97, 39.93, 31.00, 26.27, 29.88),
    'bays29': (21.31, 83.66, 98.53, 49.52, 25.38, 38.46, 28.51, 31.29, 36.54),
    'fri26': (38.54, 103.00, 90.48, 48.94, 36.56, 40.37, 43.99, 29.35, 30.52),
    'gr17': (10.28, 58.98, 50.67, 21.60, 10.44, 16.41, 16.38, 14.92, 44.84),
    'gr21': (36.58, 108.83, 78.71, 39.59, 25.12, 32.40, 32.09, 32.95, 10.68),
    'gr24': (17.22, 71.25, 80.35, 40.79, 19.00, 26.14, 21.84, 30.82, 23.11)}


def wilson(k, n, z=1.959963984540054):
    if n == 0:
        return (0.0, 1.0)
    p = k / n; d = 1 + z * z / n
    c = (p + z * z / (2 * n)) / d; h = z * math.sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / d
    return (max(0.0, c - h), min(1.0, c + h))


def load(prob, inst, method, S):
    if prob == 'mcp':
        m = method.replace(' ', '_')
        if method in ENGINE or inst in SUBSET:
            p = GSET / 'results' / 'budget' / f'G{inst}' / f'{m}_S{S}.json'
        else:
            p = HERE / 'results' / 'mcp_ext' / f'G{inst}' / f'{m}_S{S}.json'
    else:
        tag = f'G{inst}' if prob == 'gpp' else inst
        p = HERE / 'results' / prob / tag / f"{method.replace(' ', '_')}_S{S}.json"
    return json.loads(p.read_text()) if p.exists() else None


def per_run(prob, inst, d):
    """(quality, feasible, success, opt, value) arrays of the 256 final runs."""
    if prob == 'mcp':
        cuts = np.array(d['cuts']); b = BKV[f'G{inst}']
        ok = np.ones(len(cuts), bool) if d.get('finite') is None else np.array(d['finite'], bool)
        q = np.where(ok, cuts / b['BKV'], 0.0)
        return q, ok, (cuts >= b['target']) & ok, (cuts >= b['BKV']) & ok, cuts
    vals = np.array(d['values']); feas = np.array(d['feasible_runs'], bool)
    ok = np.ones(len(vals), bool) if d.get('finite') is None else np.array(d['finite'], bool)
    feas &= ok
    if prob == 'gpp':
        R = REF[f'G{inst}']['R']; tgt = REF[f'G{inst}']['target']
    else:
        R = TSPOPT[inst]; tgt = int(math.floor(1.01 * R))
    q = np.where(feas, R / np.maximum(vals, 1), 0.0)
    return q, feas, feas & (vals <= tgt), feas & (vals <= R), vals


def best_of_20(prob, q, feas):
    """Median over the 12 disjoint blocks of 20 runs (runs 0-239) of the best normalized quality in the block."""
    blocks = [q[20 * k:20 * k + 20] for k in range(12)]
    fb = [f[20 * k:20 * k + 20] for k, f in enumerate([feas] * 12)]
    best = [float(b.max()) for b in blocks]
    nofeas = sum(1 for k in range(12) if not fb[k].any())
    return float(np.median(best)), nofeas


def mcs99(S, p):
    return math.inf if p <= 0 else (float(S) if p >= 1 else S * math.log(0.01) / math.log1p(-p))


def summarize():
    out = {}
    for prob, insts in (('mcp', MCP_INST), ('gpp', GPP_INST), ('tsp', TSP_INST)):
        P = out[prob] = {}
        for inst in insts:
            key = f'G{inst}' if prob != 'tsp' else inst
            P[key] = {}
            for m in METHODS:
                rec = dict(by_S={})
                for S in S_LIST[prob]:
                    d = load(prob, inst, m, S)
                    if d is None or d.get('final') is None:
                        rec['by_S'][S] = None; continue
                    q, feas, succ, opt, vals = per_run(prob, inst, d)
                    n = len(q); k = int(succ.sum())
                    b20, nofeas = best_of_20(prob, q, feas)
                    e = dict(mean_quality=float(q.mean()), p_feasible=float(feas.mean()),
                             p_feasible_w95=wilson(int(feas.sum()), n),
                             mean_quality_feasible=float(q[feas].mean()) if feas.any() else None,
                             p_target=k / n, p_target_w95=wilson(k, n), k_opt=int(opt.sum()), best_of_20=b20,
                             best_of_20_blocks_without_feasible=nofeas, mcs99=mcs99(S, k / n),
                             selected=d.get('selected'), selected_cfg=d.get('selected_cfg'),
                             pilot_size=len(d.get('pilot', [])))
                    if prob == 'tsp' and feas.any():
                        e['arpd_feasible'] = float((100.0 * (vals[feas] - TSPOPT[inst]) / TSPOPT[inst]).mean())
                    if prob == 'mcp':
                        e['best_value'] = int(vals.max())
                    else:
                        e['best_value'] = int(vals[feas].min()) if feas.any() else None
                    if d.get('flips') is not None:
                        e['mean_flips'] = float(np.mean(d['flips']))
                    rec['by_S'][S] = e
                fin = [(v['mcs99'], S) for S, v in rec['by_S'].items() if v is not None]
                rec['mcs99_min'] = min(fin)[0] if fin else math.inf
                rec['S_mcs99_min'] = min(fin)[1] if fin and math.isfinite(min(fin)[0]) else None
                P[key][m] = rec
    return out


def fmt(x, nd=3):
    return '–' if x is None else (f'{x:.{nd}f}')


def fmt_mcs(x):
    if x is None or not math.isfinite(x):
        return '∞'
    return f'{x:,.0f}'


def geo(xs):
    xs = [x for x in xs if x is not None and math.isfinite(x) and x > 0]
    return math.exp(sum(math.log(x) for x in xs) / len(xs)) if xs else None


def tables(out):
    L = []
    for prob in ('mcp', 'gpp', 'tsp'):
        P = out[prob]; S = S_REAIM[prob]; insts = list(P)
        L.append(f'\n## {prob.upper()} (S = {S} for quality columns; MCS99 = min over S)\n')
        L.append('| Method | mean quality | feasible | quality if feasible | best-of-20 | P(target) > 0 | instances solved | '
                 'MCS99 geo-mean (solved) | runs at optimum/BKV |')
        L.append('|---|---|---|---|---|---|---|---|---|')
        for m in METHODS:
            es = [P[i][m]['by_S'].get(S) for i in insts]
            es = [e for e in es if e]
            mq = np.mean([e['mean_quality'] for e in es]); fe = np.mean([e['p_feasible'] for e in es])
            qf = [e['mean_quality_feasible'] for e in es if e['mean_quality_feasible'] is not None]
            b20 = np.mean([e['best_of_20'] for e in es])
            pt = sum(1 for e in es if e['p_target'] > 0)
            solved = sum(1 for i in insts if math.isfinite(P[i][m]['mcs99_min']))
            gm = geo([P[i][m]['mcs99_min'] for i in insts])
            kopt = sum(e['k_opt'] for e in es)
            L.append(f"| {m} | {mq:.4f} | {fe:.3f} | {fmt(np.mean(qf) if qf else None, 4)} | {b20:.4f} | {pt}/{len(insts)} | "
                     f"{solved}/{len(insts)} | {fmt_mcs(gm)} | {kopt} |")
        # per instance: quality (feasibility) at S
        L.append(f'\nPer instance, S = {S}: mean quality (feasible fraction) [P(target)]\n')
        L.append('| Inst | ' + ' | '.join(SHORT[m] for m in METHODS) + ' |')
        L.append('|---|' + '---|' * len(METHODS))
        for i in insts:
            row = []
            for m in METHODS:
                e = P[i][m]['by_S'].get(S)
                row.append('–' if not e else f"{e['mean_quality']:.3f} ({e['p_feasible']:.2f}) [{e['p_target']:.2f}]")
            L.append(f'| {i} | ' + ' | '.join(row) + ' |')
        L.append('\nPer instance: MCS99 (steps to the 1% target, min over S)\n')
        L.append('| Inst | ' + ' | '.join(SHORT[m] for m in METHODS) + ' |')
        L.append('|---|' + '---|' * len(METHODS))
        for i in insts:
            L.append(f'| {i} | ' + ' | '.join(fmt_mcs(P[i][m]['mcs99_min']) for m in METHODS) + ' |')
        if prob == 'tsp':
            L.append('\nTSP ARPD over feasible runs (%), S = 8192, with ReAIM Table VI (ASA, Neal, best of its 9 columns)\n')
            L.append('| Inst | ' + ' | '.join(SHORT[m] for m in METHODS) + ' | ReAIM ASA | Neal | ReAIM best |')
            L.append('|---|' + '---|' * (len(METHODS) + 3))
            for i in insts:
                row = [fmt(P[i][m]['by_S'][S].get('arpd_feasible'), 1) if P[i][m]['by_S'].get(S) else '–' for m in METHODS]
                t6 = REAIM_T6[i]
                L.append(f'| {i} | ' + ' | '.join(row) + f' | {t6[6]:.1f} | {t6[7]:.1f} | {min(t6):.1f} |')
        if prob == 'gpp':
            L.append('\nGPP best balanced cut (all 5 x 256 finals), with the reference R and ReAIM Table V minimum\n')
            L.append('| Inst | R | ReAIM T.V min | ' + ' | '.join(SHORT[m] for m in METHODS) + ' |')
            L.append('|---|---|---|' + '---|' * len(METHODS))
            for i in insts:
                row = []
                for m in METHODS:
                    vs = [P[i][m]['by_S'][s]['best_value'] for s in S_LIST['gpp'] if P[i][m]['by_S'].get(s) and P[i][m]['by_S'][s]['best_value'] is not None]
                    row.append(str(min(vs)) if vs else '–')
                L.append(f"| {i} | {REF[i]['R']} | {REF[i]['reaim_table_v_min']} | " + ' | '.join(row) + ' |')
    return '\n'.join(L)


def hypotheses(out):
    H = {}
    G = out['gpp']; T = out['tsp']; M = out['mcp']
    sa = {i: G[i]['SA']['by_S'][4096]['mean_quality'] for i in G}
    H['R1'] = {m: sum(1 for i in G if G[i][m]['by_S'][4096]['mean_quality'] < sa[i]) for m in ENGINE}
    H['R1_pass'] = all(v >= 8 for v in H['R1'].values())
    H['R2'] = {m: sum(1 for i in G if G[i][m]['by_S'][4096]['p_feasible'] > 0.5) for m in ('Onsager-kT', 'Onsager-online')}
    H['R2_pass'] = all(v <= 4 for v in H['R2'].values())
    best = {i: max(METHODS, key=lambda m: T[i][m]['by_S'][8192]['mean_quality']) for i in T}
    H['R3'] = dict(best_method=best, bSB_best_count=sum(1 for i in T if best[i] == 'bSB'))
    H['R3_pass'] = H['R3']['bSB_best_count'] >= 4
    r4 = {}
    for o in ('Onsager-kT', 'Onsager-online'):
        r4[o] = sum(1 for i in T if all(abs(T[i][o]['by_S'][8192]['mean_quality'] - T[i][b]['by_S'][8192]['mean_quality']) <= 0.05
                                        for b in ('SCA', 'TEC')))
    H['R4'] = r4; H['R4_pass'] = all(v >= 4 for v in r4.values())
    r5 = sum(1 for i in M if all(M[i][x]['mcs99_min'] < min(M[i][e]['mcs99_min'] for e in ENGINE) for x in ('SA', 'dSB')))
    H['R5'] = r5; H['R5_pass'] = r5 >= 15
    return H


def edges(out):
    """Grid-edge counts of the selected configurations (GPP, TSP): per method and axis, selections at the min / max."""
    sys.path.insert(0, str(HERE))
    import grids
    rep = {}
    for prob in ('gpp', 'tsp'):
        for m in METHODS:
            cnt = {}
            for i, recs in out[prob].items():
                for S, e in recs[m]['by_S'].items():
                    if not e or not e['selected_cfg']:
                        continue
                    c = e['selected_cfg']; rel = c.get('rel', {k: c[k] for k in ('dt', 'xi') if k in c})
                    for ax, v in rel.items():
                        if ax in ('kappa_factor', 'lam_factor'):
                            continue
                        cnt.setdefault(ax, []).append(v if not isinstance(v, list) else tuple(v))
            r = {}
            for ax, vs in cnt.items():
                try:
                    lo, hi = min(vs), max(vs)
                    allv = sorted(set(vs))
                except TypeError:
                    continue
                gvals = None
                r[ax] = dict(n=len(vs), distinct=len(allv), min_sel=lo, max_sel=hi)
            rep[f'{prob}/{m}'] = r
    return rep


def main():
    out = summarize()
    H = hypotheses(out)
    txt = tables(out)
    (HERE / 'results_summary.json').write_text(json.dumps(dict(summary=out, hypotheses=H), indent=1, default=str) + '\n')
    (HERE / 'analysis.txt').write_text(txt + '\n\n## Hypotheses\n\n' + json.dumps(H, indent=1, default=str) + '\n')
    print(txt)
    print(json.dumps(H, indent=1, default=str))


if __name__ == '__main__':
    main()
