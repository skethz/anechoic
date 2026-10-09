"""Addendum 4 validation summary: each baseline's reproduction of a result its paper reports, with the pre-declared
criteria of PROTOCOL_ADDENDUM4.md. Writes validation_a4_summary.json and validation_a4_summary.md."""
import json
import math
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
V = HERE / 'validation_a4'
REAIM_T4_NEAL = {'G1': 11624, 'G2': 11620, 'G6': 2177, 'G7': 2006, 'G10': 1996, 'G11': 560, 'G12': 554, 'G13': 576,
                 'G14': 3056, 'G19': 893, 'G20': 940}
REAIM_T4_ASA = {'G1': 11624, 'G2': 11615, 'G6': 2168, 'G7': 2006, 'G10': 1991, 'G11': 554, 'G12': 550, 'G13': 576,
                'G14': 3046, 'G19': 900, 'G20': 939}
REAIM_T5_ASA = {'G1': 7645, 'G2': 7611, 'G3': 7588, 'G4': 7635, 'G5': 7664, 'G14': 1115, 'G15': 1126, 'G16': 1085, 'G17': 1075}
REAIM_T6_ASA = {'bayg29': 31.00, 'bays29': 28.51, 'fri26': 43.99, 'gr17': 16.38, 'gr21': 32.09, 'gr24': 21.84}
TSPOPT = {'gr17': 2085, 'gr21': 2707, 'gr24': 1272, 'fri26': 937, 'bayg29': 1610, 'bays29': 2020}
APC_T3 = {'G22': -6545.5, 'G30': -6638.9, 'G32': -2726.3, 'G35': -3430.1}
APC_HOPT = {'G22': -6728.0, 'G30': -6826.0, 'G32': -2798.0, 'G35': -3596.0}
APC_FIG3 = {'G22': 98.3, 'G30': 97.5, 'G32': 98.05, 'G35': None}       # None: off-scale, worse than 95 %
SB21 = {('bSB', 100): 32860, ('bSB', 1000): 33215, ('bSB', 10000): 33210,
        ('dSB', 100): 32375, ('dSB', 1000): 33110, ('dSB', 10000): 33245}


def L(name):
    return json.loads((V / f'{name}.json').read_text())


def best_of_20(v, feas, maximize):
    v = np.asarray(v, float); f = np.asarray(feas, bool); out = []
    for k in range(12):
        vv = v[20 * k:20 * k + 20][f[20 * k:20 * k + 20]]
        if len(vv):
            out.append(vv.max() if maximize else vv.min())
    return float(np.median(out)) if out else None


def main():
    R = {}
    # V1a
    v1a = {}
    for inst in ('G1', 'G14', 'G22', 'K2000'):
        d = L(f'sa_vs_neal_{inst}'); a = np.array(d['ours'], float); b = np.array(d['neal'], float)
        se = math.sqrt(a.var(ddof=1) / len(a) + b.var(ddof=1) / len(b))
        v1a[inst] = dict(ours_mean=a.mean(), neal_mean=b.mean(), diff=a.mean() - b.mean(), se=se,
                         passed=bool(abs(a.mean() - b.mean()) <= 3 * se), beta_ours=d['ours_beta_range'], beta_neal=d['neal_beta_range'])
    R['V1a_SA_vs_real_Neal'] = dict(rows=v1a, passed=all(r['passed'] for r in v1a.values()))
    # V1b
    v1b = {}
    for inst, paper in REAIM_T4_NEAL.items():
        d = L(f'neal_t4_{inst}')
        bn = best_of_20(d['neal'], np.ones(len(d['neal'])), True); bo = best_of_20(d['ours'], np.ones(len(d['ours'])), True)
        v1b[inst] = dict(paper=paper, real_neal=bn, ours=bo, dev_real=(bn - paper) / paper, dev_ours=(bo - paper) / paper,
                         passed=bool(abs(bn - paper) <= 0.005 * paper and abs(bo - paper) <= 0.005 * paper))
    R['V1b_Neal_vs_ReAIM_TableIV'] = dict(rows=v1b, passed=all(r['passed'] for r in v1b.values()))
    # V2
    d = L('statica_long'); c = np.array(d['cuts'], float); se100 = c.std(ddof=1) / 10
    p = float((c >= 33000).mean())
    R['V2_STATICA_long'] = dict(mean=c.mean(), paper=33073, se100=se100, z=(c.mean() - 33073) / se100, p=p, paper_p=0.77,
                                passed=bool(abs(c.mean() - 33073) <= 3 * se100 and 0.675 <= p <= 0.848), runs=len(c))
    d = L('statica_short'); c = np.array(d['cuts'], float); se100 = c.std(ddof=1) / 10
    R['V2_STATICA_short'] = dict(mean=c.mean(), paper=32750, se100=se100, z=(c.mean() - 32750) / se100,
                                 p=float((c >= 33000).mean()), paper_p=0.07, note='known RNG-circuit effect (statica_reproduction_20261003)')
    # V3
    v3 = {}
    for jv in (0.0, 30.0, 6.0, -6.0):
        tr = []; fp = []
        for k in range(8):
            d = L(f'tec_jv{jv:g}_c{k}'); tr.append(np.array(d['cut_trace_mean'])); fp += d['first_pass']
        m = np.mean(tr, 0)
        cross = int(np.argmax(m >= 31670)) + 1 if (m >= 31670).any() else None
        fpv = np.array([x if x is not None else np.inf for x in fp], float)
        v3[jv] = dict(cycles_mean_curve=cross, median_first_pass=float(np.median(fpv)), final_mean_cut=float(m[-1]),
                      runs_reaching=int(np.isfinite(fpv).sum()), runs=len(fpv))
    t0 = v3[0.0]['cycles_mean_curve']; t30 = v3[30.0]['cycles_mean_curve']
    ok = (t0 is not None and t30 is not None and abs(t0 - 1392) <= 0.15 * 1392 and abs(t30 - 740) <= 0.15 * 740
          and abs(t0 / t30 - 1392 / 740) <= 0.2 * (1392 / 740))
    sync = {jv: float(np.mean(L(f'tec_sync_jv{jv:g}')['final_cuts'])) for jv in (0.0, 30.0)}
    R['V3_TEC_Fig3b'] = dict(rows={str(k): v for k, v in v3.items()}, paper={'0': 1392, '30': 740, 'ratio': 1392 / 740},
                             ratio=(t0 / t30) if (t0 and t30) else None, passed=bool(ok),
                             fig3a_order_ok=bool(v3[6.0]['cycles_mean_curve'] and v3[-6.0]['cycles_mean_curve'] and
                                                 v3[6.0]['cycles_mean_curve'] < t0 < v3[-6.0]['cycles_mean_curve']),
                             synchronous_misreading_final_mean_cut=sync)
    # V4
    v4a = {}
    for inst, tgt in APC_T3.items():
        e = np.array(L(f'apc_t3_{inst}')['energy_best'], float); se = e.std(ddof=1) / math.sqrt(len(e))
        v4a[inst] = dict(mean=e.mean(), paper=tgt, se=se, z=(e.mean() - tgt) / se, passed=bool(abs(e.mean() - tgt) <= 3 * se))
    R['V4a_APC_paper_SCA_Table3'] = dict(rows=v4a, passed=all(r['passed'] for r in v4a.values()))
    v4b = {}
    for inst, tgt in APC_FIG3.items():
        e = np.array(L(f'apc_fig3_{inst}')['energy_best'], float); pct = 100 * e.mean() / APC_HOPT[inst]
        v4b[inst] = dict(pct=pct, paper=tgt, passed=bool(abs(pct - tgt) <= 0.3) if tgt is not None else bool(pct < 95.0))
    R['V4b_APC_Fig3_rq045'] = dict(rows=v4b, passed=all(r['passed'] for r in v4b.values()))
    # V5
    v5 = {}
    for inst, paper in REAIM_T4_ASA.items():
        d = L(f'reaim_t4_{inst}'); b = best_of_20(d['values'], d['feasible'], True)
        v5[f'mcp/{inst}'] = dict(ours=b, paper=paper, dev=(b - paper) / paper)
    for inst, paper in REAIM_T5_ASA.items():
        d = L(f'reaim_t5_{inst}'); b = best_of_20(d['values'], d['feasible'], False)
        v5[f'gpp/{inst}'] = dict(ours=b, paper=paper, dev=(b - paper) / paper if b else None)
    for inst, paper in REAIM_T6_ASA.items():
        d = L(f'reaim_t6_{inst}'); v = np.array(d['values'], float); f = np.array(d['feasible'], bool); a = []
        for k in range(12):
            vv = v[20 * k:20 * k + 20][f[20 * k:20 * k + 20]]
            if len(vv):
                a.append(float((100 * (vv - TSPOPT[inst]) / TSPOPT[inst]).mean()))
        v5[f'tsp/{inst}'] = dict(ours_arpd=float(np.median(a)) if a else None, paper_arpd=paper, valid=float(f.mean()))
    R['V5_ReAIM_tables'] = dict(rows=v5, note='no pass/fail: k set, ITER and noise not specified by the paper')
    # V6
    d = L('asb_k2000_fig2b'); c = np.array(d['cuts'], float); fin = np.array(d['finite'], bool); cf = c[fin]
    se100 = cf.std(ddof=1) / 10
    R['V6_aSB_Fig2B'] = dict(mean=cf.mean(), paper=32768, se100=se100, z=(cf.mean() - 32768) / se100, finite=float(fin.mean()),
                             passed=bool(fin.all() and abs(cf.mean() - 32768) <= 3 * se100))
    # V7
    v7 = {}
    for (fam, S), tgt in SB21.items():
        cs = []
        for k in range({100: 1, 1000: 4, 10000: 16}[S]):
            cs += L(f'sb21_{fam}_S{S}_c{k}')['cuts']
        c = np.array(cs, float); se = c.std(ddof=1) / math.sqrt(len(c))
        v7[f'{fam}/{S}'] = dict(mean=c.mean(), read=tgt, se=se, n=len(c), passed=bool(abs(c.mean() - tgt) <= 40 + 3 * se),
                                max=float(c.max()))
    R['V7_SB_Goto2021_Fig2A'] = dict(rows=v7, passed=all(r['passed'] for r in v7.values()))
    (HERE / 'validation_a4_summary.json').write_text(json.dumps(R, indent=1, default=float) + '\n')
    print(json.dumps({k: v.get('passed') for k, v in R.items()}, indent=1))
    return R


if __name__ == '__main__':
    main()
