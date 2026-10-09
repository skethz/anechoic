#!/usr/bin/env python3
"""Pre-registered analysis of PROTOCOL_HW_MB.md (multi-bit image). Pure Python (no numpy), deterministic bootstrap.
Usage: analyze_hw_mb.py <hwmb_result_dir> <v6.4 hwv6 result dir> [--model a,b,c]
  verify/   exact gate: every trial exact_reference == PASS, all engines covered
  rounds/   K2000 regression X5, O4: TTS99 (primary over E-engine rounds, secondary from per-trial p) and R1 (same cut and
            flips per trial id as the v6.4 run) plus the per-launch cycle ratio against the v6.4 run
  gset/     G-set cohorts: TTS99 at the configuration's target, cycle-model check
Writes <dir>/mb_summary.json and prints a table. TTS conventions as analyze_hw_v6.py (round 0, which loads the couplings,
is excluded from the round time)."""
import json
import math
import random
import sys
from collections import defaultdict
from pathlib import Path


def load(p):
    return [json.loads(line) for line in Path(p).read_text().splitlines() if line.strip()]


def r99(p):
    if p >= 1: return 1.0
    if p <= 0: return math.inf
    return math.log(0.01) / math.log1p(-p)


def wilson(k, n, z=1.96):
    if n == 0: return (0.0, 1.0)
    ph = k / n; d = 1 + z * z / n; c = ph + z * z / (2 * n); h = z * math.sqrt(ph * (1 - ph) / n + z * z / (4 * n * n))
    return ((c - h) / d, (c + h) / d)


def cohort(rows, mhz, model, boot=2000, seed=20261007):
    E = rows[0]['engines']
    ok = [bool(r['success']) and bool(r['scores_agree']) for r in rows]
    p = sum(ok) / len(ok)
    rounds = defaultdict(list)
    for r, o in zip(rows, ok):
        rounds[r['round']].append((o, r['round_max_cycles'], r))
    rs = sorted(rounds)
    succ = [any(o for o, _, _ in rounds[x]) for x in rs]
    tR = [rounds[x][0][1] / (mhz * 1e3) for x in rs if x > 0]
    Pr = sum(succ) / len(succ); t = sum(tR) / len(tR)
    pred = 1 - (1 - p) ** E; se = math.sqrt(max(pred * (1 - pred), 1e-12) / len(rs))
    ind = t * (math.log(0.01) / (E * math.log1p(-p))) if 0 < p < 1 else (t if p >= 1 else math.inf)
    rnd = random.Random(seed); bt = []
    nr = len(rs)
    for _ in range(boot):
        idx = [rnd.randrange(nr) for _ in range(nr)]
        prb = sum(succ[i] for i in idx) / nr
        tb = sum((tR[i - 1] if i > 0 else t) for i in idx) / nr
        bt.append(tb if prb >= 1 else (math.inf if prb <= 0 else tb * r99(prb)))
    bt.sort()
    S = rows[0]['steps']
    # cycle model: launch cycles of one-trial launches (rounds > 0) against a*flips + b*S + c
    lc = [(r['launch_cycles'], model[0] * r['flips'] + model[1] * S + model[2]) for r in rows if r['round'] > 0 and r['launch_trials'] == 1]
    ratio = sum(a for a, _ in lc) / sum(b for _, b in lc) if lc else None
    k = sum(ok); lo, hi = wilson(k, len(ok))
    return dict(trials=len(rows), engines=E, hits=k, p=p, p_wilson95=[lo, hi], rounds=nr, P_round=Pr, P_round_pred=pred,
                M2_pass=bool(abs(Pr - pred) <= 1.96 * se + 1.0 / nr), t_round_ms=t, tts99_ms=t * r99(Pr),
                tts99_ci95=[bt[int(0.025 * boot)], bt[int(0.975 * boot) - 1]], tts99_ind_ms=ind,
                tts99_whole_ms=t * math.ceil(r99(Pr)) if Pr > 0 else math.inf, mean_flips=sum(r['flips'] for r in rows) / len(rows),
                mean_launch_cycles=sum(a for a, _ in lc) / len(lc) if lc else None, model_ratio=ratio,
                integrity_failures=sum(not r['scores_agree'] for r in rows), target=rows[0].get('target'),
                matrix=rows[0].get('matrix'), corr_scale=rows[0].get('corr_scale'), steps=S)


def main():
    d = Path(sys.argv[1]); old = Path(sys.argv[2])
    model = (0.1424, 18.09, 886.0)
    if '--model' in sys.argv: model = tuple(float(x) for x in sys.argv[sys.argv.index('--model') + 1].split(','))
    man = json.loads((d / 'image_manifest.json').read_text()); mhz = man['clock']['actual_mhz']
    out = dict(clock_mhz=mhz, logic_uuid=man['logic_uuid'], cycle_model=model, verify={}, rounds={}, gset={}, R1={})
    for f in sorted((d / 'verify').glob('*.jsonl')):
        rows = load(f)
        out['verify'][f.stem] = dict(trials=len(rows), exact=sum(r['exact_reference'] == 'PASS' for r in rows),
                                     engines=sorted({r['engine'] for r in rows}), matrix=rows[0].get('matrix'))
    out['exact_gate_pass'] = bool(out['verify']) and all(v['exact'] == v['trials'] for v in out['verify'].values())
    # K2000 exact-gate sets W1-W4 use the PROTOCOL_HW_V6 trial ids: same cut and flips as the v6.4 run's verify sets
    out['R1_verify'] = {}
    for f in sorted((d / 'verify').glob('W*.jsonl')):
        ob = old / 'verify' / f.name
        if not ob.exists(): continue
        a = {r['trial_id']: r for r in load(f)}; b = {r['trial_id']: r for r in load(ob)}
        common = sorted(set(a) & set(b))
        out['R1_verify'][f.stem] = dict(common_trials=len(common),
                                        mismatches=sum(a[i]['cut'] != b[i]['cut'] or a[i]['flips'] != b[i]['flips'] for i in common))
    for f in sorted((d / 'rounds').glob('*.jsonl')):
        rows = load(f); out['rounds'][f.stem] = cohort(rows, mhz, model)
        ob = old / 'rounds' / f.name
        if ob.exists():
            a = {r['trial_id']: r for r in rows}; b = {r['trial_id']: r for r in load(ob)}
            common = sorted(set(a) & set(b))
            bad = [i for i in common if a[i]['cut'] != b[i]['cut'] or a[i]['flips'] != b[i]['flips']]
            la = [a[i]['launch_cycles'] for i in common if a[i]['round'] > 0]; lb = [b[i]['launch_cycles'] for i in common if b[i]['round'] > 0]
            same_cyc = sum(a[i]['launch_cycles'] == b[i]['launch_cycles'] for i in common if a[i]['round'] > 0)
            out['R1'][f.stem] = dict(common_trials=len(common), mismatches=len(bad), first_mismatches=bad[:10],
                                     launch_cycle_ratio=(sum(la) / sum(lb)) if lb else None, identical_launch_cycles=same_cyc,
                                     launches_compared=len(la), v64_tts99_ms=cohort(load(ob), json.loads((old / 'image_manifest.json').read_text())['clock']['actual_mhz'], model)['tts99_ms'])
    out['R1_pass'] = bool(out['R1']) and all(v['mismatches'] == 0 and v['common_trials'] > 0 for v in out['R1'].values())
    for f in sorted((d / 'gset').glob('*.jsonl')):
        out['gset'][f.stem] = cohort(load(f), mhz, model)
    if '--predictions' in sys.argv:   # frozen PROTOCOL_HW_MB verdicts
        P = json.loads(Path(sys.argv[sys.argv.index('--predictions') + 1]).read_text()); V = {}
        V['E1'] = out['exact_gate_pass']
        V['R1'] = out['R1_pass'] and all(v['mismatches'] == 0 for v in out['R1_verify'].values()) and len(out['R1_verify']) == 4
        V['C1'] = bool(out['R1']) and all(v['launch_cycle_ratio'] is not None and abs(v['launch_cycle_ratio'] - 1) <= P['C1_tol'] for v in out['R1'].values())
        V['T1'] = all(k in out['rounds'] and abs(out['rounds'][k]['tts99_ms'] / (c['v64_tts99_ms'] * 250.0 / mhz) - 1) <= c['tol'] for k, c in P['k2000'].items())
        V['M2'] = all(v['M2_pass'] for grp in ('rounds', 'gset') for v in out[grp].values())
        for k, c in P['gset'].items():
            if k not in out['gset']: V[k] = 'missing'; continue
            g = out['gset'][k]; n1, n2 = g['trials'], c['validation_trials']
            pb = (g['hits'] + c['validation_hits']) / (n1 + n2)
            g['G1_pass'] = bool(abs(g['p'] - c['p']) <= 1.96 * math.sqrt(pb * (1 - pb) * (1 / n1 + 1 / n2)) + 1 / n1)
            F = g['mean_flips']; S = g['steps']
            tRm = (1.03 * (max(22.6 * S, 0.1424 * F + 18.09 * S) + 886) + 600) / (mhz * 1e3)
            g['G2_model_t_round_ms'] = tRm; g['G2_ratio'] = g['t_round_ms'] / tRm; g['G2_pass'] = bool(abs(g['G2_ratio'] - 1) <= P['G2_tol'])
            lo, hi = c['G3_band_ms']; g['G3_band_ms'] = [lo, hi]; g['G3_pass'] = bool(lo <= g['tts99_ms'] <= hi)
        for name in ('G1', 'G2', 'G3'):
            V[name] = all(out['gset'][k].get(name + '_pass') for k in P['gset'] if k in out['gset']) and all(k in out['gset'] for k in P['gset'])
        out['verdicts'] = V
        print('verdicts', json.dumps(V))
    (d / 'mb_summary.json').write_text(json.dumps(out, indent=2) + '\n')
    print(f"clock {mhz:.3f} MHz uuid {out['logic_uuid']}")
    for k, v in out['verify'].items(): print(f"verify {k}: exact {v['exact']}/{v['trials']} engines {len(v['engines'])} {v['matrix']}")
    print('exact gate', 'PASS' if out['exact_gate_pass'] else 'FAIL')
    for k, v in out['R1'].items(): print(f"R1 {k}: {v['mismatches']} mismatches of {v['common_trials']}; launch cycle ratio {v['launch_cycle_ratio']}; identical launch cycles {v['identical_launch_cycles']}/{v['launches_compared']}; v6.4 TTS99 {v['v64_tts99_ms']:.4f} ms")
    print('key            p      [95%]            P_round(pred)    t_round ms  TTS99 ms [95%]                 ind ms   model')
    for grp in ('rounds', 'gset'):
        for k, v in out[grp].items():
            print(f"{k:14s} {v['p']:.3f} [{v['p_wilson95'][0]:.3f},{v['p_wilson95'][1]:.3f}] {v['P_round']:.3f}({v['P_round_pred']:.3f}){'ok' if v['M2_pass'] else 'X '} "
                  f"{v['t_round_ms']:.4f}     {v['tts99_ms']:.4f} [{v['tts99_ci95'][0]:.4f},{v['tts99_ci95'][1]:.4f}]  {v['tts99_ind_ms']:.4f}  {v['model_ratio']}")


if __name__ == '__main__':
    main()
