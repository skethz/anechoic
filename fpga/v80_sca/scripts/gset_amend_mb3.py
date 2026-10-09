#!/usr/bin/env python3
"""PROTOCOL_HW_MB Amendment 3 (revision-3 image, runtime n): copy of gset_amend_mb.py (Amendment 1) for ALL instances of the
study (N = 800, 1000 and 2000; the engine runs n = N active spins with the semantics of sca_ref_bias.hpp). Keys, trial ids
and host arguments of the N = 2000 cohorts are those of Amendment 1, so they can be compared trial by trial with the
revision-1 image. 6156 trials per cohort (513 rounds of 12 engines; Amendment 1 ran 2052). Exact-gate trial ids:
950000 + 100 * (index of the instance in the sorted list of all instances).
(Amendment 1 text follows.)
PROTOCOL_HW_MB Amendment 1: turn the software G-set study's selections (research/gset_20261007/selected_configs.json)
into the frozen board configuration list. Only N = 2000 instances (the engine's N), only rules the engine implements
(plain SCA, TEC, Onsager-kT = TEC-T, Onsager-online), estimator tts_E12 (12 engines), both grid variants (original =
pre-registered, extended_amendment1); identical configurations of the two variants are run once.
Mapping to host flags (the study's model: z = s*field + q, stay clip(z/4T + 1/2), T geometric T0 -> tfin):
  plain   -> --q q --t0 T0 --t1 tfin --steps S
  tec     -> + --tec-jv jv                      (field = h + jv s(t-1); host kconst = -jv)
  tecT    -> + --tecT kappa --ramp              (field = h - kappa T(t-1) ramp s(t-1); corr scale 1)
  onsager -> + --lambda lam --ramp              (lam already includes the d/N factor -> --corr-scale 1)
Usage: gset_amend_mb.py <selected_configs.json> <out_prefix>   writes <out_prefix>_gset.txt, _verify.txt, _predictions.json"""
import json
import math
import sys

RULES = [('SCA', 'plain'), ('TEC', 'tec'), ('Onsager-kT', 'tecT'), ('Onsager-online', 'ons')]
VARS = [('original', 'o'), ('extended_amendment1', 'x')]


def args(cfg):
    f = cfg['family']
    a = f"--q {cfg['q']!r} --t0 {cfg['T0']!r} --t1 {cfg['tfin']!r} --steps {int(cfg['S'])}"
    if f == 'plain':
        a += ' --lambda 0'
    elif f == 'tec':
        a += f" --tec-jv {cfg['jv']!r}"
    elif f == 'tecT':
        assert cfg.get('ramp', False)
        a += f" --tecT {cfg['kappa']!r} --ramp"
    elif f == 'onsager':
        assert cfg.get('ramp', False)
        a += f" --lambda {cfg['lam']!r} --ramp"
    else:
        raise ValueError(f)
    return a + ' --corr-scale 1'


def main():
    src, pre = sys.argv[1], sys.argv[2]
    d = json.load(open(src))
    lines, verify, pred = [], [], {}
    insts = sorted(d['instances'], key=lambda k: int(k[1:]))   # A3: every instance (N = 800, 1000, 2000)
    for vi, inst in enumerate(insts):
        I = d['instances'][inst]
        for rule, short in RULES:
            seen = {}
            for var, vs in VARS:
                e = I['rules'].get(rule, {}).get(var, {}).get('tts_E12')
                if not e: continue
                a = args(e['cfg'])
                if a in seen:   # identical configuration in both variants: one board cohort, both labels
                    pred[seen[a]]['variants'].append(var); continue
                key = f"{inst}_{short}_{vs}"
                seen[a] = key
                n = 1024; k = e['p'] * n
                pred[key] = dict(instance=inst, N=I['N'], rule=rule, variants=[var], target=I['target'], BKV=I['BKV'], S=int(e['cfg']['S']),
                                 study_p=e['p'], study_n=n if abs(k - round(k)) < 1e-6 else None, study_hits=round(k),
                                 study_tts_E12_ms=e['tts_ms'], args=a)
                lines.append(f"{key} {inst} {I['target']} 6156 {a}")   # A3: 513 rounds of 12 (trial ids 0..2051 = Amendment 1's)
        # exact-gate set for this instance: its Onsager-online (else first available) E12 configuration, 12 trials
        for rule in ('Onsager-online', 'Onsager-kT', 'TEC', 'SCA'):
            e = I['rules'].get(rule, {}).get('original', {}).get('tts_E12')
            if e:
                verify.append(f"V{inst} {inst} {I['target']} 12 {950000 + 100 * vi} {args(e['cfg'])}"); break
    open(pre + '_gset.txt', 'w').write('# PROTOCOL_HW_MB Amendment 3 (frozen): key instance target trials host-arguments; seed 20261004, trial ids 0..trials-1,\n'
                                       '# one trial per launch, tables resident; configurations = research/gset_20261007/selected_configs.json tts_E12\n'
                                       + '\n'.join(lines) + '\n')
    open(pre + '_verify.txt', 'w').write('# PROTOCOL_HW_MB Amendment 3 exact gate (frozen): key instance target trials offset host-arguments (verify, 1 per launch)\n'
                                         + '\n'.join(verify) + '\n')
    json.dump(pred, open(pre + '_predictions.json', 'w'), indent=1)
    print(len(lines), 'cohorts,', len(verify), 'verify sets,', len(insts), 'instances')


if __name__ == '__main__':
    main()
