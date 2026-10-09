#!/usr/bin/env python3
"""Pre-registered analysis of PROTOCOL_HW_V6.md (copy of analyze_hw_multi.py with the v6 cycle model, the
pre-registered secondary estimator and the state-of-the-art comparison).
Usage: analyze_hw_v6.py <hwv6_result_dir> <cycle_model.json> <m3_ms> [sota_ms] [m6_ms]
(PROTOCOL_HW_V6: m3_ms = 0.105*300/f, sota_ms = 0.26, m6_ms = 0.05*300/f)"""
import json
import math
import sys
from collections import defaultdict
from pathlib import Path

import numpy as np

KEYS = ['P1', 'P2', 'P3', 'P4', 'T1', 'T2', 'T3', 'O1', 'O2', 'O3', 'O4']


def load(p):
    return [json.loads(line) for line in Path(p).read_text().splitlines() if line.strip()]


def r99(p):
    if p >= 1: return 1.0
    if p <= 0: return math.inf
    return math.log(0.01) / math.log1p(-p)


def main():
    d = Path(sys.argv[1])
    model = json.loads(Path(sys.argv[2]).read_text())   # cycles per trial = a*flips + b*S + c (RTL simulation fit)
    FIT = (model['a'], model['b'], model['c'])
    m3_ms = float(sys.argv[3])
    sota_ms = float(sys.argv[4]) if len(sys.argv) > 4 else 0.26   # Toshiba bSBM TTT99 (Goto 2021, measured FPGA)
    man = json.loads((d / 'image_manifest.json').read_text()); mhz = man['clock']['actual_mhz']
    S = {'clock_mhz': mhz, 'logic_uuid': man['logic_uuid'], 'verify': {}, 'rounds': {}}
    for f in sorted((d / 'verify').glob('*.jsonl')):
        rows = load(f)
        S['verify'][f.stem] = dict(trials=len(rows), exact=sum(r['exact_reference'] == 'PASS' for r in rows),
                                   engines=sorted({r['engine'] for r in rows}))
    # 2a concurrency: per-engine cycles per trial vs the v5 fit
    c = d / 'concurrency' / 'O1.jsonl'
    if c.exists():
        rows = load(c); by = defaultdict(list)
        for r in rows:
            if not r['load_J']:
                by[(r['engine'], r['round'])].append(r)
        ratios = []
        for (e, rd), g in by.items():
            pred = sum(FIT[0] * r['flips'] + FIT[1] * r['steps'] + FIT[2] for r in g)
            ratios.append(g[0]['launch_cycles'] / pred)
        S['concurrency'] = dict(launches=len(ratios), ratio_mean=float(np.mean(ratios)), ratio_min=float(np.min(ratios)),
                                ratio_max=float(np.max(ratios)), M1_pass=bool(abs(np.mean(ratios) - 1) <= 0.03),
                                p=float(np.mean([r['success'] and r['scores_agree'] for r in rows])))
    rng = np.random.default_rng(20261005)
    boots = {}
    for k in KEYS:
        f = d / 'rounds' / f'{k}.jsonl'
        if not f.exists(): continue
        rows = load(f); E = rows[0]['engines']
        ok = [r['success'] and r['scores_agree'] for r in rows]
        p = float(np.mean(ok)); n = len(rows)
        rounds = defaultdict(list)
        for r, o in zip(rows, ok):
            rounds[r['round']].append((o, r['round_max_cycles']))
        rs = sorted(rounds)
        succ = np.array([any(o for o, _ in rounds[x]) for x in rs])
        tR = np.array([rounds[x][0][1] for x in rs if x > 0]) / (mhz * 1e3)   # ms, round 0 excluded (coupling load)
        Pr = float(succ.mean()); t_mean = float(tR.mean())
        pred_Pr = 1 - (1 - p) ** E
        se = math.sqrt(max(pred_Pr * (1 - pred_Pr), 1e-12) / len(rs))
        tts = t_mean * r99(Pr); tts_whole = t_mean * math.ceil(r99(Pr)) if Pr > 0 else math.inf
        # secondary (pre-registered): P_round from the per-trial p, valid only where M2 holds
        tts_ind = t_mean * (math.log(0.01) / (E * math.log1p(-p))) if 0 < p < 1 else (t_mean if p >= 1 else math.inf)
        S['rounds'][k] = dict(trials=n, engines=E, p=p, rounds=len(rs), P_round=Pr, P_round_pred=pred_Pr,
                              M2_pass=bool(abs(Pr - pred_Pr) <= 1.96 * se + 1.0 / len(rs)),
                              t_round_ms=t_mean, t_round_max_ms=float(tR.max()), tts99_ms=tts, tts99_whole_ms=tts_whole, tts99_ind_ms=tts_ind,
                              integrity_failures=sum(not r['scores_agree'] for r in rows))
        B = 10000; nr = len(rs)
        idx = rng.integers(0, nr, size=(B, nr))
        Prb = succ[idx].mean(1)
        tb = np.concatenate([[t_mean], tR])[np.minimum(idx, len(tR))].mean(1)  # resample round times alongside
        with np.errstate(divide='ignore'):
            boots[k] = np.where(Prb >= 1, tb, np.where(Prb <= 0, np.inf, tb * np.log(0.01) / np.log1p(-np.minimum(Prb, 1 - 1e-12))))
    pt = {k: v['tts99_ms'] for k, v in S['rounds'].items()}
    if all(k in pt for k in KEYS):
        best = {fam: min((k for k in KEYS if k[0] == fam), key=lambda k: pt[k]) for fam in 'PTO'}
        comp = {}
        for fam in 'PT':
            fam_min = np.min(np.stack([boots[k] for k in KEYS if k[0] == fam]), axis=0)
            ratio = np.sort(fam_min / np.min(np.stack([boots[k] for k in KEYS if k[0] == 'O']), axis=0))
            comp[fam] = dict(best_key=best[fam], best_tts_ms=pt[best[fam]], speedup=pt[best[fam]] / pt[best['O']],
                             ci95=[float(ratio[249]), float(ratio[9749])])
        S['primary'] = dict(best_onsager=best['O'], best_onsager_tts_ms=pt[best['O']], vs=comp,
                            M3_pass=bool(pt[best['O']] <= m3_ms), M3_threshold_ms=m3_ms,
                            M4_pass=bool(comp['P']['speedup'] >= 2 and comp['T']['speedup'] >= 2),
                            M5_sota_ms=sota_ms, M5_ratio=sota_ms / pt[best['O']], M5_pass=bool(pt[best['O']] <= sota_ms / 2),
                            best_ind_onsager=min((k for k in KEYS if k[0] == 'O'), key=lambda k: S['rounds'][k]['tts99_ind_ms']),
                            best_ind_tts_ms=min(S['rounds'][k]['tts99_ind_ms'] for k in KEYS if k[0] == 'O'),
                            all_M2=bool(all(S['rounds'][k]['M2_pass'] for k in KEYS)))
        ind = {k: S['rounds'][k]['tts99_ind_ms'] for k in KEYS}
        bo = min(('O1', 'O2', 'O3', 'O4'), key=lambda k: ind[k]); bp = min(('P1', 'P2', 'P3', 'P4'), key=lambda k: ind[k])
        bt = min(('T1', 'T2', 'T3'), key=lambda k: ind[k])
        m6_ms = float(sys.argv[5]) if len(sys.argv) > 5 else 0.05
        S['primary']['M6'] = dict(best_onsager=bo, best_ms=ind[bo], vs_plain=ind[bp] / ind[bo], vs_tec=ind[bt] / ind[bo],
                                  threshold_ms=m6_ms,
                                  M6_pass=bool(ind[bo] <= m6_ms and ind[bp] / ind[bo] >= 2.5 and ind[bt] / ind[bo] >= 2.5))
    # Amendment 1 extension (re-selected schedules; compared only with re-selected baselines)
    ext = ['X1', 'X2', 'X3', 'X4', 'B1', 'B2', 'B3', 'B4'] + (['X5'] if (d / 'rounds' / 'X5.jsonl').exists() else [])
    if all((d / 'rounds' / f'{k}.jsonl').exists() for k in ext[:8]):
        for k in ext:
            rows = load(d / 'rounds' / f'{k}.jsonl'); E = rows[0]['engines']
            ok = [r['success'] and r['scores_agree'] for r in rows]; p = float(np.mean(ok))
            rounds = defaultdict(list)
            for r, o in zip(rows, ok): rounds[r['round']].append((o, r['round_max_cycles']))
            rs = sorted(rounds); succ = np.array([any(o for o, _ in rounds[x]) for x in rs])
            tR = np.array([rounds[x][0][1] for x in rs if x > 0]) / (mhz * 1e3); t = float(tR.mean()); Pr = float(succ.mean())
            pred = 1 - (1 - p) ** E; se = math.sqrt(max(pred * (1 - pred), 1e-12) / len(rs))
            S['rounds'][k] = dict(trials=len(rows), p=p, P_round=Pr, P_round_pred=pred, t_round_ms=t, tts99_ms=t * r99(Pr),
                                  tts99_ind_ms=t * (math.log(0.01) / (E * math.log1p(-p))) if 0 < p < 1 else (t if p >= 1 else math.inf),
                                  M2_pass=bool(abs(Pr - pred) <= 1.96 * se + 1.0 / len(rs)))
        R = S['rounds']
        prim_x = min(R['X1']['tts99_ms'], R['X2']['tts99_ms']); prim_b = min(R['B1']['tts99_ms'], R['B3']['tts99_ms'])
        sec_x = min(R['X3']['tts99_ind_ms'], R['X4']['tts99_ind_ms']); sec_b = min(R['B2']['tts99_ind_ms'], R['B4']['tts99_ind_ms'])
        S['amendment1'] = dict(A1_best_primary_ms=prim_x, A1_pass=bool(prim_x <= 0.08 * 300 / mhz),
                               A2_best_secondary_ms=sec_x, A2_pass=bool(sec_x <= 0.04 * 300 / mhz),
                               A3_primary_ratio=prim_b / prim_x, A3_secondary_ratio=sec_b / sec_x,
                               A3_pass=bool(prim_b / prim_x >= 2 and sec_b / sec_x >= 2.5))
        print('amendment1', json.dumps(S['amendment1']))
        if 'X5' in R and R['X5'].get('t_round_ms') and rows and rows[0]['engines'] == 12:
            # Amendment 2 (12 engines at 275 MHz): thresholds fixed for that image
            px = min(R['X1']['tts99_ms'], R['X2']['tts99_ms'], R['X5']['tts99_ms'])
            S['amendment2'] = dict(A1_best_primary_ms=px, A1_pass=bool(px <= 0.085), A2_best_secondary_ms=sec_x,
                                   A2_pass=bool(sec_x <= 0.035), A3_primary_ratio=prim_b / px, A3_secondary_ratio=sec_b / sec_x,
                                   A3_pass=bool(prim_b / px >= 2 and sec_b / sec_x >= 2.5))
            print('amendment2', json.dumps(S['amendment2']))
    (d / 'hwm_summary.json').write_text(json.dumps(S, indent=2) + '\n')
    print(f"clock {mhz:.3f} MHz uuid {S['logic_uuid']}")
    for k, v in S['verify'].items():
        print(f"verify {k}: exact {v['exact']}/{v['trials']} engines {v['engines']}")
    if 'concurrency' in S: print('concurrency', S['concurrency'])
    print('key  p      P_round (pred)   t_round ms  TTS99 ms  whole     ind')
    for k, v in S['rounds'].items():
        # (print fix after the 9-engine run: the Amendment-1 rows carry no whole-round value)
        print(f"{k:3s} {v['p']:.3f}  {v['P_round']:.3f} ({v['P_round_pred']:.3f}) {'ok' if v['M2_pass'] else 'X '}  "
              f"{v['t_round_ms']:.4f}    {v['tts99_ms']:.4f}   {v.get('tts99_whole_ms', float('nan')):.4f}   {v['tts99_ind_ms']:.4f}")
    if 'primary' in S: print('primary', json.dumps(S['primary']))


if __name__ == '__main__':
    main()
