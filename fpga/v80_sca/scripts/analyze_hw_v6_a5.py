#!/usr/bin/env python3
"""Pre-registered analysis of PROTOCOL_HW_V6 Amendment 5 (STATICA's published K2000 points on the 12 x v6.4 image).
Estimators exactly as analyze_hw_v6.py: primary TTS99 = t_R * ln(0.01) / ln(1 - P_round) with P_round the fraction of 12-engine
rounds with a success and t_R the mean device round time over rounds after round 0; secondary = t_R * ln(0.01) / (E ln(1 - p))
with p the per-trial success probability. Success = cut >= 33000 and device/host scores agree.
Verdicts: A5-C (control) X5 has the same cut and flips per trial id as the original 12 x v6.4 X5 run (2052 of 2052);
A5-M (reported) SP_long / SP_short p within 1.96 sqrt(pbar (1 - pbar)(2/2052)) + 0.02 of the reference-model p
(research/fairness_audit_20261008/k2000_pub: 0.829 and 0.138).
Ratios: primary SP / O1 (O1 = X5, this session and the original run), secondary SP / O5 (O5 = X4, original run).
Usage: analyze_hw_v6_a5.py <hwv6_..._a5 dir> <original hwv6_board_v64_e12_250mhz dir>"""
import json
import math
import sys
from collections import defaultdict
from pathlib import Path

MODEL_P = {'SP_long': (0.8289473684210527, 2052), 'SP_short': (0.13791423001949318, 2052)}
MODEL_MS = {'SP_long': (0.30846031112307204, 0.06703878853819048), 'SP_short': (0.24856772907220842, 0.24856772907220848)}


def load(p):
    return [json.loads(line) for line in Path(p).read_text().splitlines() if line.strip()]


def r99(p):
    if p >= 1: return 1.0
    if p <= 0: return math.inf
    return math.log(0.01) / math.log1p(-p)


def est(rows, mhz):
    E = rows[0]['engines']
    ok = [r['success'] and r['scores_agree'] for r in rows]; p = sum(ok) / len(ok)
    rounds = defaultdict(list)
    for r, o in zip(rows, ok): rounds[r['round']].append((o, r['round_max_cycles']))
    rs = sorted(rounds); succ = [any(o for o, _ in rounds[x]) for x in rs]
    tR = [rounds[x][0][1] / (mhz * 1e3) for x in rs if x > 0]; t = sum(tR) / len(tR); Pr = sum(succ) / len(succ)
    pred = 1 - (1 - p) ** E; se = math.sqrt(max(pred * (1 - pred), 1e-12) / len(rs))
    return dict(trials=len(rows), engines=E, hits=sum(ok), p=p, rounds=len(rs), P_round=Pr, P_round_pred=pred,
                M2_pass=bool(abs(Pr - pred) <= 1.96 * se + 1.0 / len(rs)), t_round_ms=t, tts99_ms=t * r99(Pr),
                tts99_ind_ms=t * (math.log(0.01) / (E * math.log1p(-p))) if 0 < p < 1 else (t if p >= 1 else math.inf),
                mean_flips=sum(r['flips'] for r in rows) / len(rows), integrity_failures=sum(not r['scores_agree'] for r in rows))


def main():
    d, o = Path(sys.argv[1]), Path(sys.argv[2])
    mhz = json.loads((d / 'image_manifest.json').read_text())['clock']['actual_mhz']
    mhz0 = json.loads((o / 'image_manifest.json').read_text())['clock']['actual_mhz']
    S = dict(clock_mhz=mhz, rounds={}, original={})
    for k in ('X5', 'SP_long', 'SP_short'):
        f = d / 'rounds' / f'{k}.jsonl'
        if f.exists(): S['rounds'][k] = est(load(f), mhz)
    for k in ('X5', 'X4', 'P3'):
        S['original'][k] = est(load(o / 'rounds' / f'{k}.jsonl'), mhz0)
    a = {r['trial_id']: r for r in load(d / 'rounds' / 'X5.jsonl')}; b = {r['trial_id']: r for r in load(o / 'rounds' / 'X5.jsonl')}
    common = sorted(set(a) & set(b)); bad = [i for i in common if a[i]['cut'] != b[i]['cut'] or a[i]['flips'] != b[i]['flips']]
    S['A5_C'] = dict(common=len(common), mismatches=len(bad), pass_=bool(len(common) == 2052 and not bad))
    S['A5_M'] = {}
    for k, (pm, nm) in MODEL_P.items():
        if k not in S['rounds']: continue
        c = S['rounds'][k]; nb = c['trials']; pb = (c['hits'] + pm * nm) / (nb + nm)
        S['A5_M'][k] = dict(p_board=c['p'], p_model=pm, within=bool(abs(c['p'] - pm) <= 1.96 * math.sqrt(pb * (1 - pb) * (1 / nb + 1 / nm)) + 0.02),
                            model_primary_ms=MODEL_MS[k][0], model_secondary_ms=MODEL_MS[k][1])
    R, O = S['rounds'], S['original']
    S['ratios'] = {}
    for k in ('SP_long', 'SP_short'):
        if k not in R: continue
        S['ratios'][k] = dict(primary_SP_over_O1_session=R[k]['tts99_ms'] / R['X5']['tts99_ms'],
                              primary_SP_over_O1_original=R[k]['tts99_ms'] / O['X5']['tts99_ms'],
                              secondary_SP_over_O5_original=R[k]['tts99_ind_ms'] / O['X4']['tts99_ind_ms'],
                              primary_O1_session_over_SP=R['X5']['tts99_ms'] / R[k]['tts99_ms'],
                              secondary_O5_original_over_SP=O['X4']['tts99_ind_ms'] / R[k]['tts99_ind_ms'])
    (d / 'a5_summary.json').write_text(json.dumps(S, indent=1) + '\n')
    print('A5-C X5 control identical to the original run:', S['A5_C'])
    print('A5-M p against the reference model:', json.dumps(S['A5_M']))
    print('key       p      P_round (pred)    t_round ms  TTS99 primary ms  secondary ms  mean flips')
    for k, c in list(R.items()) + [(f'orig {k}', v) for k, v in O.items()]:
        print(f"{k:9s} {c['p']:.4f} {c['P_round']:.4f} ({c['P_round_pred']:.4f}) {'ok' if c['M2_pass'] else 'X '} {c['t_round_ms']:.4f}      "
              f"{c['tts99_ms']:.4f}            {c['tts99_ind_ms']:.4f}       {c['mean_flips']:.0f}")
    for k, v in S['ratios'].items(): print('ratios', k, {x: round(y, 3) for x, y in v.items()})


if __name__ == '__main__':
    main()
