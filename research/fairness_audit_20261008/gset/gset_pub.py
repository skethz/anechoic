"""PROTOCOL_GSET_PUB.md: Table 8a baselines (plain SCA, TEC) at the settings their papers describe, engine model only.
No board action: the selected configurations are pre-registered for a later board run and the script stops.

Rules (absolute units; G-set and K2000 both have |J| = 1, and neither paper gives a G-set setting or a scaling recipe):
  SCA-pub   STATICA (Yamamoto et al., JSSC 2021): clipped Eq. (7), q = 4, T_fin = 5, geometric T; STATICA's K2000
            grid S in {360, 560, ..., 1560} x T_init in {30, 40, 50}, chosen by TTS (here the 12-engine primary TTS of
            Table 8a, v6.4 cycle model), on 256-run pilots.
  TEC-pub   TEC's published J_v = +30 J (Du et al., arXiv:2608.21753, K2000 optimum) added to SCA-pub's rule
            (field h + J_v s(t-1)), same grid and selection. TEC's own dynamics (sequential Glauber SA) cannot run
            on the engine.
  SCA-apc   context only: the SCA recipe for G-set max-cut of the APC-SCA paper (Okonogi et al., IEICE 2023, by the
            STATICA group): q = lambda/2 (lambda = largest eigenvalue of -J), exponential T 10 -> 0.1, S = 1000;
            the engine's clipped rule replaces the paper's logistic one.
Confirmation: the selected configuration, 1,024 fresh runs. Seeds: pilots SeedSequence([20261008, 501, g, rule, k]),
confirmations [20261008, 502, g, rule], SCA-apc [20261008, 503, g].
Usage: python3 gset_pub.py run [workers] | summarize"""
import json
import math
import os
import sys
import time
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
AUDIT = Path(os.environ.get('AUDIT_ROOT', HERE.parents[2]))
for sub in ('research/theory_ideas_20261003', 'research/statica_reproduction_20261003', 'research/ablation_20261005',
            'research/reaim_reproduction_20261003', 'research/algorithm_compare_20261007', 'research/gset_20261007'):
    if str(AUDIT / sub) not in sys.path:
        sys.path.append(str(AUDIT / sub))
import engine as E  # noqa: E402   (research/gset_20261007/engine.py, unedited; bit-identical to abl.run)

GS = AUDIT / 'research/gset_20261007'
OUT = HERE / 'results_gset_pub'
SEED = 20261008
S_GRID = (360, 560, 760, 960, 1160, 1360, 1560)
T_INIT = (30.0, 40.0, 50.0)
Q, T_FIN, JV = 4.0, 5.0, 30.0
B_PILOT, B_CONF = 256, 1024
ALL = tuple(range(1, 48)) + (51, 52, 53, 54)
RULES = {'SCA-pub': 0, 'TEC-pub': 1}


def grid(rule):
    out = []
    for S in S_GRID:
        for T0 in T_INIT:
            c = dict(family='plain', q=Q, T0=T0, S=S)
            if rule == 'TEC-pub':
                c = dict(family='tec', jv=JV, q=Q, T0=T0, S=S)
            out.append(c)
    return out


def bkv():
    return json.loads((GS / 'bkv.json').read_text())['instances']


def tts12(cuts, flips, S, target):
    p = float((np.asarray(cuts) >= target).mean())
    tr = E.expected_max(E.t_trial_ms(np.asarray(flips), S), E.E_ENGINES); P = 1 - (1 - p) ** E.E_ENGINES
    prim = tr if P >= 0.995 else (math.inf if p <= 0 else tr * E.r99(P))
    sec = math.inf if p <= 0 else (tr if p >= 1 else tr * math.log(0.01) / (E.E_ENGINES * math.log1p(-p)))
    return dict(p=p, t_round12_ms=tr, P_round12=P, tts12_ms=prim, tts12_secondary_ms=sec,
                tts1_ms=E.t_trial_ms(float(np.mean(flips)), S) * E.r99(p))


_G = {}


def graph(g):
    if g not in _G:
        _G[g] = E.Graph(f'G{g}')
    return _G[g]


def simulate(g, cfg, B, seed, tfin):
    G = graph(g)
    rng = np.random.default_rng(np.random.SeedSequence(seed))
    s0 = rng.choice(np.array([-1.0, 1.0], dtype=np.float32), size=(B, G.N))
    s, flips = E.run(G, s0, cfg, rng, tfin)
    return G.cut(s), flips


def job(args):
    kind, g, rule, k, cfg, seed, B, tfin = args
    t0 = time.time()
    target = bkv()[f'G{g}']['target']
    cuts, flips = simulate(g, cfg, B, seed, tfin)
    r = dict(kind=kind, g=g, rule=rule, k=k, cfg=cfg, seed=seed, runs=B, tfin=tfin, target=target,
             mean_cut=float(cuts.mean()), max_cut=int(cuts.max()), mean_flips=float(flips.mean()),
             **tts12(cuts, flips, cfg['S'], target), sec=time.time() - t0)
    if kind != 'pilot':
        r['cuts'] = cuts.tolist(); r['flips'] = flips.tolist()
    return r


def path(r):
    return OUT / 'raw' / f"G{r['g']}_{r['rule']}_{r['kind']}_{r['k']}.json"


def select(pilots):
    """Minimum 12-engine primary TTS; if infinite for all, the highest mean cut; ties to the lower grid index."""
    return min(pilots, key=lambda r: (r['tts12_ms'], -r['mean_cut'], r['k']))


def lam_half(g):
    G = graph(g)
    return float(np.linalg.eigvalsh(-G.dense().astype(np.float64))[-1]) / 2


def run(workers):
    (OUT / 'raw').mkdir(parents=True, exist_ok=True)
    pilots = [('pilot', g, rule, k, c, [SEED, 501, g, RULES[rule], k], B_PILOT, T_FIN)
              for g in ALL for rule in RULES for k, c in enumerate(grid(rule))]
    apc = []
    for g in ALL:
        apc.append(('apc', g, 'SCA-apc', 0, dict(family='plain', q=lam_half(g), T0=10.0, S=1000), [SEED, 503, g], B_CONF, 0.1))
    todo = [a for a in pilots + apc if not path(dict(g=a[1], rule=a[2], kind=a[0], k=a[3])).exists()]
    print(len(todo), 'pilot/apc jobs', flush=True)
    t0 = time.time()
    with ProcessPoolExecutor(max_workers=workers) as ex:
        futs = [ex.submit(job, a) for a in todo]
        for n, f in enumerate(as_completed(futs), 1):
            r = f.result(); path(r).write_text(json.dumps(r) + '\n')
            if n % 100 == 0:
                print(f'{n}/{len(todo)} {time.time() - t0:.0f}s', flush=True)
        conf = []
        for g in ALL:
            for rule in RULES:
                pil = [json.loads(path(dict(g=g, rule=rule, kind='pilot', k=k)).read_text()) for k in range(len(grid(rule)))]
                s = select(pil)
                a = ('conf', g, rule, s['k'], s['cfg'], [SEED, 502, g, RULES[rule]], B_CONF, T_FIN)
                if not path(dict(g=g, rule=rule, kind='conf', k=s['k'])).exists():
                    conf.append(a)
        print(len(conf), 'confirmation jobs', flush=True)
        for r in ex.map(job, conf):
            path(r).write_text(json.dumps(r) + '\n')
    summarize()


def summarize():
    sel = json.loads((GS / 'selected_configs.json').read_text())['instances']
    rows = {}
    for g in ALL:
        row = {}
        for rule in RULES:
            pil = [json.loads(path(dict(g=g, rule=rule, kind='pilot', k=k)).read_text()) for k in range(len(grid(rule)))]
            s = select(pil)
            c = json.loads(path(dict(g=g, rule=rule, kind='conf', k=s['k'])).read_text())
            row[rule] = {x: c[x] for x in ('cfg', 'p', 'tts12_ms', 'tts12_secondary_ms', 'tts1_ms', 't_round12_ms', 'mean_cut', 'target')}
            row[rule]['pilot_p_max'] = max(r['p'] for r in pil)
        a = json.loads(path(dict(g=g, rule='SCA-apc', kind='apc', k=0)).read_text())
        row['SCA-apc'] = {x: a[x] for x in ('cfg', 'p', 'tts12_ms', 'tts12_secondary_ms', 'tts1_ms', 't_round12_ms', 'mean_cut', 'target')}
        for var in ('original', 'extended_amendment1'):
            for rule in ('SCA', 'TEC', 'Onsager-kT', 'Onsager-online'):
                v = sel[f'G{g}']['rules'][rule][var]['tts_E12']
                row[f'{rule} ({var})'] = dict(tts12_ms=v['tts_ms'], p=v['p'], S=v['S'])
        rows[g] = row
    OUT.mkdir(exist_ok=True)
    (OUT / 'per_instance.json').write_text(json.dumps(rows, indent=1) + '\n')
    table(rows)


CLASSES = {'Random, +1': list(range(1, 6)) + list(range(22, 27)) + list(range(43, 48)),
           'Random, +-1': list(range(6, 11)) + list(range(27, 32)),
           'Toroidal, +-1': [11, 12, 13, 32, 33, 34],
           'Planar-like, +1': [14, 15, 16, 17, 35, 36, 37, 38, 51, 52, 53, 54],
           'Planar-like, +-1': [18, 19, 20, 21, 39, 40, 41, 42]}


def table(rows):
    lines = ['Model 12-engine TTS99 speedup = TTS(baseline)/TTS(corrected rule); geometric mean over instances where both are '
             'finite; (wins / n); [baseline never reaches the target in the 1,024-run confirmation on k instances]']
    out = {}
    for var in ('original', 'extended_amendment1'):
        for base in ('SCA-pub', 'TEC-pub', 'SCA-apc', f'SCA ({var})', f'TEC ({var})'):
            for corr in (f'Onsager-kT ({var})', f'Onsager-online ({var})'):
                for cls, inst in list(CLASSES.items()) + [('All 51', ALL)]:
                    r = [(rows[g][base]['tts12_ms'], rows[g][corr]['tts12_ms']) for g in inst]
                    fin = [a / b for a, b in r if math.isfinite(a) and math.isfinite(b)]
                    wins = sum(1 for a, b in r if a > b)
                    never = sum(1 for a, b in r if not math.isfinite(a))
                    gm = math.exp(sum(map(math.log, fin)) / len(fin)) if fin else float('nan')
                    out[f'{base} vs {corr} | {cls}'] = dict(geomean=gm, wins=wins, n=len(inst), baseline_never=never)
                    lines.append(f'{base:16s} vs {corr:33s} {cls:17s} {gm:7.2f} ({wins}/{len(inst)}) [{never}]')
    for rule in ('SCA-pub', 'TEC-pub', 'SCA-apc'):
        never = [g for g in ALL if rows[g][rule]['p'] <= 0]
        lines.append(f'{rule}: target never reached in confirmation on {len(never)}/51: {never}')
    t = '\n'.join(lines)
    print(t)
    (OUT / 'table8a_model.txt').write_text(t + '\n')
    (OUT / 'table8a_model.json').write_text(json.dumps(out, indent=1) + '\n')


if __name__ == '__main__':
    if sys.argv[1] == 'run':
        run(int(sys.argv[2]) if len(sys.argv) > 2 else 48)
    else:
        summarize()
