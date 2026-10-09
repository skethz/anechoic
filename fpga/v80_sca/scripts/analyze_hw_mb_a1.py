#!/usr/bin/env python3
"""Pre-registered analysis of PROTOCOL_HW_MB Amendment 1 (the G-set study's 12-engine selections on the multi-bit image).
Usage: analyze_hw_mb_a1.py <hwmb_..._a1 dir> PROTOCOL_HW_MB_A1_predictions.json
Per cohort: TTS99 as analyze_hw_mb.py (cohort()); A1-P: |p_board - p_study| <= 1.96 sqrt(pbar(1-pbar)(1/n_b + 1/n_s)) + 0.02;
A1-T: board primary TTS99 / study tts_E12 (reported, median ratio predicted in [0.8, 1.5]). Writes a1_summary.json."""
import json
import math
import statistics as st
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from analyze_hw_mb import cohort, load  # noqa: E402


def main():
    d = Path(sys.argv[1]); P = json.loads(Path(sys.argv[2]).read_text())
    mhz = json.loads((d / 'image_manifest.json').read_text())['clock']['actual_mhz']
    model = (0.1424, 18.09, 886.0)
    out = dict(clock_mhz=mhz, verify={}, cohorts={})
    for f in sorted((d / 'verify').glob('*.jsonl')):
        rows = load(f)
        out['verify'][f.stem] = dict(trials=len(rows), exact=sum(r['exact_reference'] == 'PASS' for r in rows),
                                     engines=len({r['engine'] for r in rows}))
    out['A1_E_pass'] = len(out['verify']) == 21 and all(v['exact'] == v['trials'] == 12 for v in out['verify'].values())
    for key, pr in P.items():
        f = d / 'gset' / f'{key}.jsonl'
        if not f.exists(): out['cohorts'][key] = 'missing'; continue
        c = cohort(load(f), mhz, model)
        nb, ns = c['trials'], pr['study_n']; pb = (c['hits'] + pr['study_hits']) / (nb + ns)
        c['study_p'] = pr['study_p']; c['study_tts_E12_ms'] = pr['study_tts_E12_ms']
        c['A1_P_pass'] = bool(abs(c['p'] - pr['study_p']) <= 1.96 * math.sqrt(pb * (1 - pb) * (1 / nb + 1 / ns)) + 0.02)
        c['tts_ratio_vs_study'] = c['tts99_ms'] / pr['study_tts_E12_ms'] if pr['study_tts_E12_ms'] else None
        for k in ('instance', 'rule', 'variants', 'S', 'target'): c[k] = pr[k]
        out['cohorts'][key] = c
    done = [c for c in out['cohorts'].values() if isinstance(c, dict)]
    ratios = [c['tts_ratio_vs_study'] for c in done if c['tts_ratio_vs_study'] and math.isfinite(c['tts_ratio_vs_study'])]
    out['A1_P_fraction'] = sum(c['A1_P_pass'] for c in done) / max(1, len(done))
    out['A1_P_pass'] = len(done) == len(P) and out['A1_P_fraction'] >= 0.90
    out['A1_T_median_ratio'] = st.median(ratios) if ratios else None
    out['A1_T_pass'] = bool(ratios) and 0.8 <= out['A1_T_median_ratio'] <= 1.5
    # best measured 12-engine TTS per instance (any rule; both variants)
    best = {}
    for k, c in out['cohorts'].items():
        if not isinstance(c, dict): continue
        b = best.get(c['instance'])
        if b is None or c['tts99_ms'] < b[1]: best[c['instance']] = (k, c['tts99_ms'])
    out['best_per_instance'] = {i: dict(key=k, tts99_ms=t) for i, (k, t) in sorted(best.items(), key=lambda x: int(x[0][1:]))}
    (d / 'a1_summary.json').write_text(json.dumps(out, indent=1) + '\n')
    print('A1-E exact gate', out['A1_E_pass'], {k: f"{v['exact']}/{v['trials']}" for k, v in out['verify'].items()})
    print(f"A1-P fraction {out['A1_P_fraction']:.3f} pass {out['A1_P_pass']}; A1-T median TTS ratio {out['A1_T_median_ratio']} pass {out['A1_T_pass']}")
    for k, c in out['cohorts'].items():
        if isinstance(c, dict):
            print(f"{k:18s} S {c['S']:5d} p {c['p']:.3f} (study {c['study_p']:.3f}) {'ok' if c['A1_P_pass'] else 'X '} "
                  f"t_R {c['t_round_ms']:.4f} TTS99 {c['tts99_ms']:.4f} ms (study {c['study_tts_E12_ms']:.4f}) M2 {'ok' if c['M2_pass'] else 'X'}")
    for i, b in out['best_per_instance'].items(): print('best', i, b)


if __name__ == '__main__':
    main()
