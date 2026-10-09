#!/usr/bin/env python3
"""PROTOCOL_HW_MB Amendment 4 lists from the ReAIM hardware export (research/reaim_benchmarks_20261008/hw_export/
selected_configs.json, copied with hashes to data/a4_problems): functional checks of the bias path on the K = 4 + bias image.
  verify: every exported TSP instance (bias, n = 289..841, S = 8192) with its SCA schedule, and three GPP instances (dense
          K = 4, b = 0) with their SCA schedule: 6 trials each (one per engine), every trial bit-exact against run_trial_bias;
  cohorts: every exported TSP instance with its SCA schedule, 1026 trials (trial ids 0..1025, the export's seed), to compare
          the number of valid tours (permutation matrices) among trial ids 0..1023 with the export's precision_study.
Usage: a4_lists.py <selected_configs.json> <data dir on fpga-host> <out prefix>"""
import json
import sys


def args(c):
    a = f"--t0 {c['t0']!r} --t1 {c['t1']!r} --steps {int(c['S'])} --q {c['q']!r}"
    if c['family'] == 'plain': a += ' --lambda 0'
    elif c['family'] == 'tec': a += f" --tec-jv {c['jv']!r}"
    elif c['family'] == 'tecT': a += f" --tecT {c['kappa']!r}" + (' --ramp' if c['ramp'] else '')
    elif c['family'] == 'onsager': a += f" --lambda {c['lam']!r}" + (' --ramp' if c['ramp'] else '')
    else: raise ValueError(c['family'])
    return a + ' --corr-scale 1 --bias-mode 1'


def main():
    src, ddir, pre = sys.argv[1:4]
    d = json.load(open(src))['instances']
    ver, coh, pred = [], [], {}
    k = 0
    for name, v in d.items():
        e = v['K']['4']; c = e['rules']['SCA']
        if not all(c['checks'][x] for x in c['checks'] if x != 'max_T'): continue
        spec = f"jint8:{ddir}/{e['jint8']}:{ddir}/{e['bias']}"
        if v['problem'] == 'tsp' or name in ('G1', 'G14', 'G17'):
            ver.append(f"V{name} {spec} -1000000000 6 {970000 + 100 * k} --seed {c['seed']} {args(c)}"); k += 1
        if v['problem'] == 'tsp':
            key = f"{name}_plain"
            coh.append(f"{key} {spec} -1000000000 1026 --seed {c['seed']} {args(c)}")
            ps = c['precision_study']
            pred[key] = dict(problem=name, kind='tsp', n=v['N'], S=int(c['S']), seed=c['seed'], ref_trials=ps['trials'],
                             ref_p_feasible=ps['p_feasible'], ref_feasible=round(ps['p_feasible'] * ps['trials']),
                             ref_p_target=ps['p_target'], ref_mean_quality=ps['mean_quality'], args=args(c))
    open(pre + '_verify.txt', 'w').write('# PROTOCOL_HW_MB Amendment 4 exact gate (frozen): key problem target_energy trials offset host-arguments\n' + '\n'.join(ver) + '\n')
    open(pre + '_cohorts.txt', 'w').write('# PROTOCOL_HW_MB Amendment 4 cohorts (frozen): key problem target_energy trials host-arguments (trial ids 0..trials-1)\n' + '\n'.join(coh) + '\n')
    json.dump(dict(cohorts=pred), open(pre + '_predictions.json', 'w'), indent=1)
    print(len(ver), 'verify sets,', len(coh), 'cohorts')


if __name__ == '__main__':
    main()
