#!/usr/bin/env python3
"""Pre-registered analysis of PROTOCOL_HW_MB Amendment 4 (revision-3 K = 4 + bias image, 6 engines).
Usage: analyze_hw_mb_a4.py <hwmb_..._a4 dir> <v6.4 hwv6 dir> PROTOCOL_HW_MB_A4_predictions.json
Verdicts (PROTOCOL_HW_MB_A4.md):
  A4-E  every exact-gate trial bit-exact against run_trial_bias, padding bits 0, every set covers all engines
  A4-R  K2000 with bias mode off: X5 and O4 same cut and flips per trial id as the v6.4 run, W1-W4 identical to v6.4's sets
  A4-C  launch-cycle ratio against v6.4 within 1% (X5, O4)
  A4-F  per TSP cohort: the number of valid tours (the n = c^2 one-hot spins form a permutation matrix) among trial ids
        0..1023 equals the export's precision_study count (p_feasible x 1024) exactly - the export used the same golden model,
        seed and trial ids, so a bit-exact engine reproduces it trial for trial
Reported only (functional checks, no quality claims): valid-tour fraction over all 1026 trials, best and median energy,
round time. Writes a4_summary.json."""
import json
import math
import struct
import statistics as st
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from analyze_hw_mb import cohort, load  # noqa: E402


def same(a_rows, b_rows):
    a = {r['trial_id']: r for r in a_rows}; b = {r['trial_id']: r for r in b_rows}
    common = sorted(set(a) & set(b))
    bad = [i for i in common if a[i]['cut'] != b[i]['cut'] or a[i]['flips'] != b[i]['flips']]
    return common, bad, a, b


def main():
    d, v64, P = Path(sys.argv[1]), Path(sys.argv[2]), json.loads(Path(sys.argv[3]).read_text())
    mhz = json.loads((d / 'image_manifest.json').read_text())['clock']['actual_mhz']
    mhz64 = json.loads((v64 / 'image_manifest.json').read_text())['clock']['actual_mhz']
    model = (0.1424, 18.09, 886.0)
    out = dict(clock_mhz=mhz, verify={}, R={}, rounds={}, cohorts={})
    E = None
    for f in sorted((d / 'verify').glob('*.jsonl')):
        rows = load(f); E = rows[0]['engines']
        out['verify'][f.stem] = dict(trials=len(rows), exact=sum(r['exact_reference'] == 'PASS' for r in rows),
                                     engines=len({r['engine'] for r in rows}), n=rows[0].get('n'), bias_mode=rows[0].get('bias_mode'),
                                     padding_zero=all(r.get('padding_zero', False) for r in rows),
                                     best_energy=min((r['energy'] for r in rows if 'energy' in r), default=None))
    nver = 10 + sum(1 for l in (d / 'PROTOCOL_HW_MB_A4_verify.txt').read_text().splitlines() if l.strip() and not l.startswith('#'))
    out['A4_E'] = (len(out['verify']) == nver and all(v['exact'] == v['trials'] and v['engines'] == E and v['padding_zero']
                                                      for v in out['verify'].values()))
    okR = okC = True
    for k in ('W1', 'W2', 'W3', 'W4'):
        a, b = d / 'verify' / f'{k}.jsonl', v64 / 'verify' / f'{k}.jsonl'
        if not (a.exists() and b.exists()): okR = False; out['R'][k] = 'missing'; continue
        common, bad, _, _ = same(load(a), load(b)); out['R'][k] = dict(common=len(common), mismatches=len(bad))
        okR &= bool(common) and not bad
    for k in ('X5', 'O4'):
        a, b = d / 'rounds' / f'{k}.jsonl', v64 / 'rounds' / f'{k}.jsonl'
        if not (a.exists() and b.exists()): okR = okC = False; out['R'][k] = 'missing'; continue
        ra, rb = load(a), load(b)
        common, bad, A, B = same(ra, rb)
        both = [i for i in common if A[i]['round'] > 0 and B[i]['round'] > 0]   # launches without the coupling load in both runs
        la = sum(A[i]['launch_cycles'] for i in both); lb = sum(B[i]['launch_cycles'] for i in both)
        ca, cb = cohort(ra, mhz, model), cohort(rb, mhz64, model)
        out['rounds'][k] = ca
        out['R'][k] = dict(common=len(common), mismatches=len(bad), first=bad[:10], launch_cycle_ratio=la / lb if lb else None,
                           tts99_ms=ca['tts99_ms'], v64_tts99_ms=cb['tts99_ms'], engines=ca['engines'])
        okR &= len(common) == len(rb) and not bad
        okC &= out['R'][k]['launch_cycle_ratio'] is not None and abs(out['R'][k]['launch_cycle_ratio'] - 1) <= 0.01
    out['A4_R'] = okR; out['A4_C'] = okC
    okF = True
    for key, pr in P.get('cohorts', {}).items():
        f = d / 'cohorts' / f'{key}.jsonl'
        if not f.exists(): out['cohorts'][key] = 'missing'; okF = False; continue
        rows = load(f)
        raw = (d / 'cohorts' / f'{key}.spins.bin').read_bytes()
        n = rows[0]['n']; cc = int(round(math.sqrt(n)))
        valid = {}
        for k, r in enumerate(rows):   # spins.bin holds 32 words per trial, in the order of the jsonl lines
            w = struct.unpack_from('<32Q', raw, 256 * k)
            bit = lambda i: (w[i >> 6] >> (i & 63)) & 1
            ok = cc * cc == n and all(sum(bit(v * cc + j) for j in range(cc)) == 1 and sum(bit(j * cc + v) for j in range(cc)) == 1 for v in range(cc))
            valid[r['trial_id']] = ok
        c = cohort(rows, mhz, model)
        en = sorted(r['energy'] for r in rows)
        nref = sum(valid.get(i, False) for i in range(pr['ref_trials']))
        c.update(best_energy=en[0], median_energy=en[len(en) // 2], n=n, valid_tours_all=sum(valid.values()), trials_all=len(valid),
                 valid_tours_ref_ids=nref, ref_valid_tours=pr['ref_feasible'], ref_p_feasible=pr['ref_p_feasible'],
                 ref_mean_quality=pr['ref_mean_quality'], ref_p_target=pr['ref_p_target'],
                 A4_F_pass=bool(nref == pr['ref_feasible'] and all(i in valid for i in range(pr['ref_trials']))))
        okF &= c['A4_F_pass']
        out['cohorts'][key] = c
    out['A4_F'] = okF and bool(P.get('cohorts'))
    out['verdicts'] = {k: out[k] for k in ('A4_E', 'A4_R', 'A4_C', 'A4_F')}
    out['gset'] = {k: dict(tts99_ms=c['tts99_ms']) for k, c in out['cohorts'].items() if isinstance(c, dict)}   # for analyze_power_mb.py
    (d / 'a4_summary.json').write_text(json.dumps(out, indent=1) + '\n')
    print('verdicts', json.dumps(out['verdicts']))
    for k, v in out['verify'].items(): print(f"verify {k}: {v['exact']}/{v['trials']} engines {v['engines']} n {v['n']} bias {v['bias_mode']} best E {v['best_energy']}")
    for k, v in out['R'].items(): print('R', k, v)
    for k, c in out['cohorts'].items():
        if isinstance(c, dict):
            print(f"{k:16s} n {c['n']:4d} S {c['steps']:5d} valid tours ids 0..1023: {c['valid_tours_ref_ids']} (export {c['ref_valid_tours']}) "
                  f"{'ok' if c['A4_F_pass'] else 'X'}; all {c['valid_tours_all']}/{c['trials_all']}; t_R {c['t_round_ms']:.4f} ms; best E {c['best_energy']}")


if __name__ == '__main__':
    main()
