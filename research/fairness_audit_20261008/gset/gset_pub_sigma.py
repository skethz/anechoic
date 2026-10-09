"""PROTOCOL_GSET_PUB_A1.md: the published K2000 settings transferred to each G-set instance with the COP study's rule
(research/reaim_benchmarks_20261008/PROTOCOL_ADDENDUM3.md: K2000 values multiplied by alpha = sigma/sigma_K2000,
sigma = sqrt(mean_i sum_j J_ij^2), sigma_K2000 = sqrt(1999); the rule is ours, as there labelled).
  SCA-pub-sigma  STATICA long point: q = 4 alpha, T 40 alpha -> 5 alpha
  TEC-pub-sigma  the same plus TEC's J_v = 30 alpha (engine adaptation: field h + J_v s(t-1))
S from the budget set of every Table 8a rule, {250, 500, 1000, 2000, 4000}, chosen by the 12-engine primary TTS99 on
256-run pilots (ties / all infinite -> highest mean cut, then smaller S); 1,024-run confirmation.
Seeds: pilots [20261008, 511, g, rule, S index], confirmations [20261008, 512, g, rule]. Usage: run [workers] | summarize"""
import json
import math
import sys
import time
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
import gset_pub as GP  # noqa: E402  (hashed; imported unchanged)

E = GP.E
OUT = GP.HERE / 'results_gset_pub_sigma'
S_LIST = (250, 500, 1000, 2000, 4000)
RULES = {'SCA-pub-sigma': 0, 'TEC-pub-sigma': 1}
SIGMA_K = math.sqrt(1999.0)


def cfg(g, rule, S):
    G = GP.graph(g); a = G.sigma / SIGMA_K
    c = dict(family='plain', q=4.0 * a, T0=40.0 * a, S=S)
    if rule == 'TEC-pub-sigma':
        c = dict(family='tec', jv=30.0 * a, q=4.0 * a, T0=40.0 * a, S=S)
    return c, 5.0 * a


def path(kind, g, rule, k):
    return OUT / 'raw' / f'G{g}_{rule}_{kind}_{k}.json'


def run(workers):
    (OUT / 'raw').mkdir(parents=True, exist_ok=True)
    jobs = []
    for g in GP.ALL:
        for rule in RULES:
            for si, S in enumerate(S_LIST):
                c, tf = cfg(g, rule, S)
                if not path('pilot', g, rule, si).exists():
                    jobs.append(('pilot', g, rule, si, c, [GP.SEED, 511, g, RULES[rule], si], GP.B_PILOT, tf))
    print(len(jobs), 'pilot jobs', flush=True)
    with ProcessPoolExecutor(max_workers=workers) as ex:
        futs = [ex.submit(GP.job, a) for a in jobs]
        for n, f in enumerate(as_completed(futs), 1):
            r = f.result(); path('pilot', r['g'], r['rule'], r['k']).write_text(json.dumps(r) + '\n')
            if n % 100 == 0:
                print(n, flush=True)
        conf = []
        for g in GP.ALL:
            for rule in RULES:
                pil = [json.loads(path('pilot', g, rule, si).read_text()) for si in range(len(S_LIST))]
                s = min(pil, key=lambda r: (r['tts12_ms'], -r['mean_cut'], r['k']))
                if not path('conf', g, rule, s['k']).exists():
                    conf.append(('conf', g, rule, s['k'], s['cfg'], [GP.SEED, 512, g, RULES[rule]], GP.B_CONF, s['tfin']))
        print(len(conf), 'confirmations', flush=True)
        for r in ex.map(GP.job, conf):
            path('conf', r['g'], r['rule'], r['k']).write_text(json.dumps(r) + '\n')
    summarize()


def summarize():
    sel = json.loads((GP.GS / 'selected_configs.json').read_text())['instances']
    rows = {}
    for g in GP.ALL:
        row = {}
        for rule in RULES:
            pil = [json.loads(path('pilot', g, rule, si).read_text()) for si in range(len(S_LIST))]
            s = min(pil, key=lambda r: (r['tts12_ms'], -r['mean_cut'], r['k']))
            c = json.loads(path('conf', g, rule, s['k']).read_text())
            row[rule] = {x: c[x] for x in ('cfg', 'tfin', 'p', 'tts12_ms', 'tts12_secondary_ms', 'tts1_ms', 't_round12_ms', 'mean_cut', 'target')}
        for var in ('original', 'extended_amendment1'):
            for rule in ('SCA', 'TEC', 'Onsager-kT', 'Onsager-online'):
                v = sel[f'G{g}']['rules'][rule][var]['tts_E12']
                row[f'{rule} ({var})'] = dict(tts12_ms=v['tts_ms'], p=v['p'], S=v['S'])
        rows[g] = row
    (OUT / 'per_instance.json').write_text(json.dumps(rows, indent=1) + '\n')
    lines = ['Model 12-engine TTS99 speedup = TTS(baseline)/TTS(corrected); geometric mean over instances with both finite; '
             '(wins/n) [baseline never reaches target]']
    out = {}
    for var in ('original', 'extended_amendment1'):
        for base in ('SCA-pub-sigma', 'TEC-pub-sigma'):
            for corr in (f'Onsager-kT ({var})', f'Onsager-online ({var})'):
                for cls, inst in list(GP.CLASSES.items()) + [('All 51', GP.ALL)]:
                    r = [(rows[g][base]['tts12_ms'], rows[g][corr]['tts12_ms']) for g in inst]
                    fin = [a / b for a, b in r if math.isfinite(a) and math.isfinite(b)]
                    gm = math.exp(sum(map(math.log, fin)) / len(fin)) if fin else float('nan')
                    wins = sum(1 for a, b in r if a > b); never = sum(1 for a, b in r if not math.isfinite(a))
                    out[f'{base} vs {corr} | {cls}'] = dict(geomean=gm, wins=wins, n=len(inst), baseline_never=never)
                    lines.append(f'{base:14s} vs {corr:35s} {cls:17s} {gm:7.2f} ({wins}/{len(inst)}) [{never}]')
    for rule in RULES:
        never = [g for g in GP.ALL if rows[g][rule]['p'] <= 0]
        lines.append(f'{rule}: target never reached in confirmation on {len(never)}/51: {never}')
    t = '\n'.join(lines)
    print(t)
    (OUT / 'table8a_model_sigma.txt').write_text(t + '\n')
    (OUT / 'table8a_model_sigma.json').write_text(json.dumps(out, indent=1) + '\n')


if __name__ == '__main__':
    if sys.argv[1] == 'run':
        run(int(sys.argv[2]) if len(sys.argv) > 2 else 32)
    else:
        summarize()
