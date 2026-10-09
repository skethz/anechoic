#!/usr/bin/env python3
"""Pre-registered analysis of PROTOCOL_HW_MULTI.md. Usage: analyze_hw_multi.py <hwm_result_dir> [m3_ms m3_ref_mhz]
M3 threshold defaults to the original 0.20 ms x (300/f); Amendment 3 (3 engines) uses 0.25 ms x (250/f)."""
import json
import math
import sys
from collections import defaultdict
from pathlib import Path

import numpy as np

V5_FIT = (0.2754, 31.06, 753.0)   # v5 measured cycles per run (couplings resident, batched)
KEYS = ['P1', 'P2', 'P3', 'P4', 'T1', 'T2', 'T3', 'O1', 'O2', 'O3', 'O4']


def load(p):
    return [json.loads(line) for line in Path(p).read_text().splitlines() if line.strip()]


def r99(p):
    if p >= 1: return 1.0
    if p <= 0: return math.inf
    return math.log(0.01) / math.log1p(-p)


def main():
    d = Path(sys.argv[1])
    m3_ms = float(sys.argv[2]) if len(sys.argv) > 2 else 0.20
    m3_ref = float(sys.argv[3]) if len(sys.argv) > 3 else 300.0
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
            pred = sum(V5_FIT[0] * r['flips'] + V5_FIT[1] * r['steps'] + V5_FIT[2] for r in g)
            ratios.append(g[0]['launch_cycles'] / pred)
        S['concurrency'] = dict(launches=len(ratios), ratio_mean=float(np.mean(ratios)), ratio_min=float(np.min(ratios)),
                                ratio_max=float(np.max(ratios)), M1_pass=bool(abs(np.mean(ratios) - 1) <= 0.02),
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
        S['rounds'][k] = dict(trials=n, engines=E, p=p, rounds=len(rs), P_round=Pr, P_round_pred=pred_Pr,
                              M2_pass=bool(abs(Pr - pred_Pr) <= 1.96 * se + 1.0 / len(rs)),
                              t_round_ms=t_mean, t_round_max_ms=float(tR.max()), tts99_ms=tts, tts99_whole_ms=tts_whole,
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
                            M3_pass=bool(pt[best['O']] <= m3_ms * m3_ref / mhz), M3_threshold_ms=m3_ms * m3_ref / mhz,
                            M4_pass=bool(comp['P']['speedup'] >= 2 and comp['T']['speedup'] >= 2))
    (d / 'hwm_summary.json').write_text(json.dumps(S, indent=2) + '\n')
    print(f"clock {mhz:.3f} MHz uuid {S['logic_uuid']}")
    for k, v in S['verify'].items():
        print(f"verify {k}: exact {v['exact']}/{v['trials']} engines {v['engines']}")
    if 'concurrency' in S: print('concurrency', S['concurrency'])
    print('key  p      P_round (pred)   t_round ms  TTS99 ms  whole')
    for k, v in S['rounds'].items():
        print(f"{k:3s} {v['p']:.3f}  {v['P_round']:.3f} ({v['P_round_pred']:.3f}) {'ok' if v['M2_pass'] else 'X '}  "
              f"{v['t_round_ms']:.4f}    {v['tts99_ms']:.4f}   {v['tts99_whole_ms']:.4f}")
    if 'primary' in S: print('primary', json.dumps(S['primary']))


if __name__ == '__main__':
    main()
