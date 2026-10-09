#!/usr/bin/env python3
"""Pre-registered analysis of PROTOCOL_HW_MB Amendment 3 (revision-3 K = 2 image with runtime n).
Usage: analyze_hw_mb_a3.py <hwmb_..._a3 dir> <v6.4 hwv6 dir> <revision-1 Amendment-1 dir> PROTOCOL_HW_MB_A3_predictions.json
Verdicts (PROTOCOL_HW_MB_A3.md):
  A3-E  every exact-gate trial bit-exact (all sets; every set covers all engines)
  A3-R  K2000 X5 and O4: same cut and flips per trial id as the v6.4 run (2052 each) and W1-W4 identical to v6.4's sets
  A3-C  launch-cycle ratio against v6.4 within 1% (X5, O4)
  A3-I  the 100 N = 2000 cohorts: trial ids 0..2051 have the same cut and flips as Amendment 1 on the revision-1 image
  A3-P  per cohort |p_board - p_study| <= 1.96 sqrt(pbar (1 - pbar) (1/n_b + 1/n_s)) + 0.02; pass if >= 90% of cohorts
  A3-T  median over cohorts of board primary TTS99 / study tts_E12 in [0.8, 1.5]
Reported (not verdicts): best measured 12-engine TTS99 per instance (post-hoc minimum over rules and variants), geometric
means per rule and instance class. TTS estimators: analyze_hw_mb.cohort (primary over rounds, round 0 excluded).
Writes a3_summary.json."""
import json
import math
import statistics as st
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from analyze_hw_mb import cohort, load  # noqa: E402


def gclass(inst):
    i = int(inst[1:])
    if i <= 10 or 22 <= i <= 31 or 43 <= i <= 47: return 'random'
    if 11 <= i <= 13 or 32 <= i <= 34: return 'toroidal'
    return 'planar'


def same(a_rows, b_rows, ids=None):
    a = {r['trial_id']: r for r in a_rows}; b = {r['trial_id']: r for r in b_rows}
    common = sorted(set(a) & set(b)) if ids is None else [i for i in ids if i in a and i in b]
    bad = [i for i in common if a[i]['cut'] != b[i]['cut'] or a[i]['flips'] != b[i]['flips']]
    return common, bad, a, b


def gmean(xs):
    xs = [x for x in xs if x and math.isfinite(x) and x > 0]
    return math.exp(sum(map(math.log, xs)) / len(xs)) if xs else None


def main():
    d, v64, r1, P = Path(sys.argv[1]), Path(sys.argv[2]), Path(sys.argv[3]), json.loads(Path(sys.argv[4]).read_text())
    mhz = json.loads((d / 'image_manifest.json').read_text())['clock']['actual_mhz']
    mhz64 = json.loads((v64 / 'image_manifest.json').read_text())['clock']['actual_mhz']
    model = (0.1424, 18.09, 886.0)
    out = dict(clock_mhz=mhz, verify={}, R={}, I={}, cohorts={})
    E = None
    for f in sorted((d / 'verify').glob('*.jsonl')):
        rows = load(f); E = rows[0]['engines']
        out['verify'][f.stem] = dict(trials=len(rows), exact=sum(r['exact_reference'] == 'PASS' for r in rows),
                                     engines=len({r['engine'] for r in rows}), n=rows[0].get('n'),
                                     padding_zero=all(r.get('padding_zero', False) for r in rows))
    nver = 10 + sum(1 for l in (d / 'PROTOCOL_HW_MB_A3_verify.txt').read_text().splitlines() if l.strip() and not l.startswith('#'))
    out['A3_E'] = (len(out['verify']) == nver and all(v['exact'] == v['trials'] and v['engines'] == E and v['padding_zero']
                                                      for v in out['verify'].values()))
    # K2000 regression against v6.4
    okR = True; okC = True
    for k in ('W1', 'W2', 'W3', 'W4'):
        a, b = d / 'verify' / f'{k}.jsonl', v64 / 'verify' / f'{k}.jsonl'
        if not (a.exists() and b.exists()): okR = False; out['R'][k] = 'missing'; continue
        common, bad, _, _ = same(load(a), load(b)); out['R'][k] = dict(common=len(common), mismatches=len(bad))
        okR &= bool(common) and not bad
    out['rounds'] = {}
    for k in ('X5', 'O4'):
        a, b = d / 'rounds' / f'{k}.jsonl', v64 / 'rounds' / f'{k}.jsonl'
        if not (a.exists() and b.exists()): okR = okC = False; out['R'][k] = 'missing'; continue
        ra, rb = load(a), load(b)
        common, bad, A, B = same(ra, rb)
        la = sum(A[i]['launch_cycles'] for i in common if A[i]['round'] > 0); lb = sum(B[i]['launch_cycles'] for i in common if B[i]['round'] > 0)
        ca, cb = cohort(ra, mhz, model), cohort(rb, mhz64, model)
        out['rounds'][k] = ca
        out['R'][k] = dict(common=len(common), mismatches=len(bad), first=bad[:10], launch_cycle_ratio=la / lb if lb else None,
                           identical_launch_cycles=sum(A[i]['launch_cycles'] == B[i]['launch_cycles'] for i in common if A[i]['round'] > 0),
                           tts99_ms=ca['tts99_ms'], v64_tts99_ms=cb['tts99_ms'])
        okR &= len(common) == len(rb) and not bad
        okC &= out['R'][k]['launch_cycle_ratio'] is not None and abs(out['R'][k]['launch_cycle_ratio'] - 1) <= 0.01
    out['A3_R'] = okR; out['A3_C'] = okC
    # G-set cohorts
    okI = True; nI = 0
    for key, pr in P.items():
        f = d / 'gset' / f'{key}.jsonl'
        if not f.exists(): out['cohorts'][key] = 'missing'; continue
        rows = load(f)
        c = cohort(rows, mhz, model)
        nb, ns = c['trials'], pr['study_n']; pb = (c['hits'] + pr['study_hits']) / (nb + ns)
        c['study_p'] = pr['study_p']; c['study_tts_E12_ms'] = pr['study_tts_E12_ms']
        c['A3_P_pass'] = bool(abs(c['p'] - pr['study_p']) <= 1.96 * math.sqrt(pb * (1 - pb) * (1 / nb + 1 / ns)) + 0.02)
        c['tts_ratio_vs_study'] = c['tts99_ms'] / pr['study_tts_E12_ms'] if pr['study_tts_E12_ms'] else None
        for k in ('instance', 'N', 'rule', 'variants', 'S', 'target', 'BKV'): c[k] = pr[k]
        c['class'] = gclass(pr['instance'])
        c['best_cut'] = max(r['cut'] for r in rows); c['best_cut_over_BKV'] = c['best_cut'] / pr['BKV']
        if pr['N'] == 2000:
            g = r1 / 'gset' / f'{key}.jsonl'
            if not g.exists(): okI = False; c['I'] = 'missing'
            else:
                common, bad, _, _ = same(rows, load(g), ids=range(2052))
                c['I'] = dict(common=len(common), mismatches=len(bad)); okI &= len(common) == 2052 and not bad; nI += 1
        out['cohorts'][key] = c
    done = [c for c in out['cohorts'].values() if isinstance(c, dict)]
    out['A3_I'] = okI and nI == sum(1 for p in P.values() if p['N'] == 2000)
    ratios = [c['tts_ratio_vs_study'] for c in done if c['tts_ratio_vs_study'] and math.isfinite(c['tts_ratio_vs_study'])]
    out['A3_P_fraction'] = sum(c['A3_P_pass'] for c in done) / max(1, len(done))
    out['A3_P'] = len(done) == len(P) and out['A3_P_fraction'] >= 0.90
    out['A3_T_median_ratio'] = st.median(ratios) if ratios else None
    out['A3_T'] = bool(ratios) and 0.8 <= out['A3_T_median_ratio'] <= 1.5
    best = {}
    for k, c in out['cohorts'].items():
        if not isinstance(c, dict): continue
        b = best.get(c['instance'])
        if b is None or c['tts99_ms'] < b['tts99_ms']:
            best[c['instance']] = dict(key=k, rule=c['rule'], S=c['S'], N=c['N'], cls=c['class'], p=c['p'], tts99_ms=c['tts99_ms'],
                                       tts99_ci95=c['tts99_ci95'], study_tts_E12_ms=c['study_tts_E12_ms'])
    out['best_per_instance'] = dict(sorted(best.items(), key=lambda x: int(x[0][1:])))
    gm = {}
    for rule in ('SCA', 'TEC', 'Onsager-kT', 'Onsager-online'):
        for cls in ('random', 'toroidal', 'planar', 'all'):
            per = {}
            for c in done:
                if c['rule'] != rule or (cls != 'all' and c['class'] != cls) or 'original' not in c['variants']: continue
                per[c['instance']] = c['tts99_ms']
            gm[f'{rule}/{cls}'] = dict(instances=len(per), geomean_tts99_ms=gmean(per.values()))
    for cls in ('random', 'toroidal', 'planar', 'all'):
        v = [b['tts99_ms'] for i, b in best.items() if cls == 'all' or b['cls'] == cls]
        gm[f'best/{cls}'] = dict(instances=len(v), geomean_tts99_ms=gmean(v))
    out['geomeans'] = gm
    out['verdicts'] = {k: out[k] for k in ('A3_E', 'A3_R', 'A3_C', 'A3_I', 'A3_P', 'A3_T')}
    out['gset'] = {k: dict(tts99_ms=c['tts99_ms']) for k, c in out['cohorts'].items() if isinstance(c, dict)}   # for analyze_power_mb.py
    (d / 'a3_summary.json').write_text(json.dumps(out, indent=1) + '\n')
    print('verdicts', json.dumps(out['verdicts']))
    print(f"A3-P fraction {out['A3_P_fraction']:.3f}; A3-T median ratio {out['A3_T_median_ratio']}")
    for k, v in out['verify'].items(): print(f"verify {k}: {v['exact']}/{v['trials']} engines {v['engines']} n {v['n']} padding0 {v['padding_zero']}")
    for k, v in out['R'].items(): print('R', k, v)
    for k, c in out['cohorts'].items():
        if isinstance(c, dict):
            print(f"{k:18s} N {c['N']:4d} S {c['S']:5d} p {c['p']:.3f} (study {c['study_p']:.3f}) {'ok' if c['A3_P_pass'] else 'X '} "
                  f"t_R {c['t_round_ms']:.4f} TTS99 {c['tts99_ms']:.4f} ms (study {c['study_tts_E12_ms']:.4f}) M2 {'ok' if c['M2_pass'] else 'X'}"
                  + (f" I {c['I']}" if 'I' in c else ''))
    for i, b in out['best_per_instance'].items(): print('best', i, b)
    for k, v in gm.items(): print('geomean', k, v)


if __name__ == '__main__':
    main()
