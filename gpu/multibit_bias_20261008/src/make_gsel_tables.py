#!/usr/bin/env python3
"""Amendment 5 analysis: GPU-selected G-set schedules against the V80's own measured A3 results.
Inputs (on gpu-host, project root): results/GS/plan_labels.json, results/GS/plan_cohorts.tsv, results/GSR (held-out cohorts),
results/GSW (power), data/v80_a3/a3_summary.json and power_summary.json (unchanged copies of the V80 study's files),
data/gset_study_results/{conf,ext_conf} (the study's confirmation runs, for the labelled V80 modelled secondary), results/G
(the earlier identical-schedule runs). Estimators: src/analyze.py (unchanged). Writes results/GS_tables.{md,json}."""
import datetime, glob, json, math, os, statistics, sys
sys.path.insert(0, 'src')
import analyze as A

RULES = ['SCA', 'TEC', 'Onsager-kT', 'Onsager-online']
CLASS = {}
for g in list(range(1, 11)) + list(range(22, 32)) + list(range(43, 48)): CLASS[f'G{g}'] = 'random'
for g in list(range(11, 14)) + list(range(32, 35)): CLASS[f'G{g}'] = 'toroidal'
for g in list(range(14, 22)) + list(range(35, 43)) + list(range(51, 55)): CLASS[f'G{g}'] = 'planar'
CLASSES = ['random', 'toroidal', 'planar']
LN01 = math.log(0.01)

def fmt(x, d=4):
    if x is None: return '—'
    if isinstance(x, float) and math.isinf(x): return '∞'
    return f'{x:.{d}f}'

def gmean(xs):
    xs = list(xs)
    if not xs or any(x is None or x <= 0 or math.isinf(x) for x in xs): return None
    return math.exp(sum(math.log(x) for x in xs) / len(xs))

labels = json.load(open('results/GS/plan_labels.json'))
cohorts = {l.split('\t')[0]: l.rstrip('\n').split('\t') for l in open('results/GS/plan_cohorts.tsv')}
insts = sorted({x['instance'] for x in labels}, key=lambda g: int(g[1:]))
est = {}
for key in cohorts:
    pre = f'results/GSR/{key}'
    if os.path.exists(pre + '.summary.json'):
        e = A.estimate(pre); s = json.load(open(pre + '.summary.json'))
        e.update(variant=f"R{s['r']}G{s['groups']}", mismatches=s['mismatches'])
        tr = [json.loads(l) for l in open(pre + '.trials.jsonl')]
        e['mean_cut'] = sum(t['cut'] for t in tr) / len(tr); e['best_cut'] = max(t['cut'] for t in tr)
        est[key] = e
out = dict(n_cohorts=len(cohorts), n_measured=len(est), mismatches=sum(e['mismatches'] for e in est.values()))
md = []
# ---------------- verification ----------------
if os.path.exists('results/GSV/summary.json'):
    v = json.load(open('results/GSV/summary.json')); out['verification'] = v
    md.append(f"## Amendment 5: GPU-selected G-set schedules\n\nVerification: {v['exact']}/{v['trials']} traced trials bit-exact against run_trial_bias in {v['runs']} runs "
              f"(one per unique selected cohort); device/host mismatches {v['device_host_mismatches']}; final-field mismatches {v['field_mismatches']}.\n")
md.append(f"Held-out cohorts: {len(est)}/{len(cohorts)} measured, 128 batches each; device/host cut mismatches {out['mismatches']}.\n")
# ---------------- GPU per label ----------------
G = {}
for x in labels:
    co = x['cohort']; e = est.get(co)
    if not e: continue
    G[(x['instance'], x['rule'], x['grid'], x['estimator'])] = dict(cohort=co, S=x['S'], B=x['B'], p=e['p'], P_batch=e['P_batch'], t_ms=e['t_batch_ms'],
        primary=e['tts_primary_ms'], primary_ci=e['tts_primary_ci_ms'], secondary=e['tts_secondary_ms'], secondary_ci=e['tts_secondary_ci_ms'],
        pred_primary=x['pred_primary_ms'], pred_secondary=x['pred_secondary_ms'], variant=e['variant'], mean_cut=e['mean_cut'], best_cut=e['best_cut'],
        BKV=x['BKV'], v80sec=x.get('v80_secondary_selection'))
out['GPU'] = {'|'.join(k): v for k, v in G.items()}
# ---------------- V80 A3 (unchanged copy) ----------------
a3 = json.load(open('data/v80_a3/a3_summary.json'))['cohorts']
V = {}
for key, c in a3.items():
    for var in c['variants']:
        grid = 'original' if var == 'original' else 'extended'
        V[(c['instance'], c['rule'], grid)] = dict(key=key, S=c['S'], p=c['p'], P_round=c['P_round'], t_round_ms=c['t_round_ms'], primary=c['tts99_ms'],
                                                   primary_ci=c['tts99_ci95'], secondary=c['tts99_ind_ms'])
out['V80'] = {'|'.join(k): v for k, v in V.items()}
# V80 secondary at its own secondary-optimal selection (modelled: study float-model p, V80 cycle model t_round12)
def v80_model_secondary(g, r, grid):
    sel = G.get((g, r, grid, 'secondary'), {}).get('v80sec')
    if not sel: return None, None
    S = sel['S']; cd = 'conf' if grid == 'original' else 'ext_conf'
    best = (sel['k'], sel['n'], sel['t_round12_ms'], 'final256')
    for obj in ('E1', 'E12'):
        f = f'data/gset_study_results/{cd}/{g}/{r}_{obj}.json'
        if os.path.exists(f):
            d = json.load(open(f))
            if d['S'] == S and d['final'].get('t_round12_ms') is not None:
                best = (d['final']['k_target'], d['final']['runs'], d['final']['t_round12_ms'], f'conf{obj}1024'); break
    k, n, tr, src = best
    if k == 0: return math.inf, src
    p = (n - 0.5) / n if k == n else k / n
    return tr * LN01 / (12 * math.log(1 - p)), src

# ---------------- comparison tables ----------------
cmp = {}
for grid in ('original', 'extended'):
    md.append(f"\n### {grid.capitalize()} grid: GPU (own selection) against V80 (own selection, A3 board)\n")
    md.append('| Rule | GPU primary (ms) | V80 primary (ms) | V80/GPU | GPU faster (of 51) | GPU secondary (ms) | V80 secondary, A3 at E12 (ms) | V80/GPU | GPU faster | V80 secondary, own selection, modelled (ms) | V80/GPU |')
    md.append('|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|')
    for r in RULES + ['best of 4 rules']:
        rows = []
        for g in insts:
            rr = RULES if r == 'best of 4 rules' else [r]
            gp = min(G[(g, x, grid, 'primary')]['primary'] for x in rr); vp = min(V[(g, x, grid)]['primary'] for x in rr)
            gs = min(G[(g, x, grid, 'secondary')]['secondary'] for x in rr); vs = min(V[(g, x, grid)]['secondary'] for x in rr)
            vm = min(v80_model_secondary(g, x, grid)[0] for x in rr)
            rows.append(dict(instance=g, cls=CLASS[g], gp=gp, vp=vp, gs=gs, vs=vs, vm=vm))
        G_p, V_p, G_s, V_s, V_m = (gmean(x[k] for x in rows) for k in ('gp', 'vp', 'gs', 'vs', 'vm'))
        wp = sum(x['gp'] < x['vp'] for x in rows); ws = sum(x['gs'] < x['vs'] for x in rows)
        cmp[f'{grid}|{r}'] = dict(rows=rows, gpu_primary=G_p, v80_primary=V_p, gpu_secondary=G_s, v80_secondary_a3=V_s, v80_secondary_model=V_m,
                                 gpu_wins_primary=wp, gpu_wins_secondary=ws, gpu_wins_secondary_vs_model=sum(x['gs'] < x['vm'] for x in rows),
                                 per_class={c: dict(n=sum(x['cls'] == c for x in rows), gpu_primary=gmean(x['gp'] for x in rows if x['cls'] == c),
                                                    v80_primary=gmean(x['vp'] for x in rows if x['cls'] == c), gpu_secondary=gmean(x['gs'] for x in rows if x['cls'] == c),
                                                    v80_secondary_a3=gmean(x['vs'] for x in rows if x['cls'] == c), v80_secondary_model=gmean(x['vm'] for x in rows if x['cls'] == c),
                                                    gpu_wins_primary=sum(x['gp'] < x['vp'] for x in rows if x['cls'] == c),
                                                    gpu_wins_secondary=sum(x['gs'] < x['vs'] for x in rows if x['cls'] == c)) for c in CLASSES})
        ratio = lambda a, b: (a / b) if (a and b) else None
        md.append(f"| {r} | {fmt(G_p)} | {fmt(V_p)} | {fmt(ratio(V_p, G_p), 2)} | {wp} | {fmt(G_s, 5)} | {fmt(V_s, 5)} | {fmt(ratio(V_s, G_s), 2)} | {ws} | {fmt(V_m, 5)} | {fmt(ratio(V_m, G_s), 2)} |")
    md.append(f"\nPer graph class ({grid} grid; geometric means in ms; V80/GPU in brackets; GPU wins of n):\n")
    md.append('| Rule | Class | n | GPU primary | V80 primary | GPU wins | GPU secondary | V80 secondary A3 | GPU wins | V80 secondary modelled |')
    md.append('|---|---|---:|---:|---:|---:|---:|---:|---:|---:|')
    for r in RULES + ['best of 4 rules']:
        for c in CLASSES:
            x = cmp[f'{grid}|{r}']['per_class'][c]
            q = lambda a, b: f" [{a / b:.2f}]" if (a and b) else ''
            md.append(f"| {r} | {c} | {x['n']} | {fmt(x['gpu_primary'])} | {fmt(x['v80_primary'])}{q(x['v80_primary'], x['gpu_primary'])} | {x['gpu_wins_primary']} | "
                      f"{fmt(x['gpu_secondary'], 5)} | {fmt(x['v80_secondary_a3'], 5)}{q(x['v80_secondary_a3'], x['gpu_secondary'])} | {x['gpu_wins_secondary']} | {fmt(x['v80_secondary_model'], 5)}{q(x['v80_secondary_model'], x['gpu_secondary'])} |")
out['comparison'] = cmp
# ---------------- earlier identical-schedule comparison (E12 cells, best of B), for reference ----------------
if os.path.exists('results/G/plan.json'):
    plan = json.load(open('results/G/plan.json')); old = {}
    for c in plan:
        if c['kind'] != 'E12': continue
        es = [A.estimate(f"results/G/{c['tag']}_b{B}") for B in (60, 120, 240) if os.path.exists(f"results/G/{c['tag']}_b{B}.summary.json")]
        old[(c['instance'], c['rule'])] = (min(e['tts_primary_ms'] for e in es), min(e['tts_secondary_ms'] for e in es))
    # Amendment 6: phase G ran the study's extended-grid (extended_amendment1) E12 configurations, so the identical-schedule
    # pairing is with the V80's extended-grid cohorts. Two B conventions: per-instance best of B, and fixed B (60 for the
    # primary, 240 for the secondary, the best B at the geometric-mean level, which is the convention of the paper's figures).
    for c in plan:
        if c['kind'] != 'E12': continue
        for B in (60, 120, 240):
            e = A.estimate(f"results/G/{c['tag']}_b{B}"); old[(c['instance'], c['rule'], B)] = (e['tts_primary_ms'], e['tts_secondary_ms'])
    md.append('\n### For reference: the earlier identical-schedule comparison (GPU on the study\'s extended-grid E12 configurations of phase G, against the V80\'s extended-grid A3 cohorts, which are the same configurations)\n')
    md.append('| Rule | GPU primary, best B per instance (ms) | GPU primary, B = 60 (ms) | V80 primary (ms) | V80/GPU (best B; B = 60) | GPU secondary, best B (ms) | GPU secondary, B = 240 (ms) | V80 secondary (ms) | V80/GPU (best B; B = 240) |')
    md.append('|---|---:|---:|---:|---|---:|---:|---:|---|')
    oldt = {}
    for r in RULES:
        gp = gmean(old[(g, r)][0] for g in insts); gs = gmean(old[(g, r)][1] for g in insts)
        gp60 = gmean(old[(g, r, 60)][0] for g in insts); gs240 = gmean(old[(g, r, 240)][1] for g in insts)
        vp = gmean(V[(g, r, 'extended')]['primary'] for g in insts); vs = gmean(V[(g, r, 'extended')]['secondary'] for g in insts)
        oldt[r] = dict(gpu_primary_bestB=gp, gpu_primary_B60=gp60, v80_primary=vp, gpu_secondary_bestB=gs, gpu_secondary_B240=gs240, v80_secondary=vs)
        md.append(f"| {r} | {fmt(gp)} | {fmt(gp60)} | {fmt(vp)} | {vp / gp:.2f}; {vp / gp60:.2f} | {fmt(gs, 5)} | {fmt(gs240, 5)} | {fmt(vs, 5)} | {vs / gs:.2f}; {vs / gs240:.2f} |")
    out['identical_schedule_reference'] = oldt
# ---------------- GPU selection details ----------------
md.append('\n### GPU selections (original grid): distribution of (S, B)\n')
from collections import Counter
for e_ in ('primary', 'secondary'):
    c = Counter((x['S'], x['B']) for x in labels if x['grid'] == 'original' and x['estimator'] == e_)
    md.append(f"- {e_}: " + ', '.join(f"S={S} B={B}: {n}" for (S, B), n in sorted(c.items())))
pr = [(G[k]['primary'], G[k]['pred_primary']) for k in G if k[3] == 'primary']
sr = [(G[k]['secondary'], G[k]['pred_secondary']) for k in G if k[3] == 'secondary']
rat = lambda xs: statistics.median(a / b for a, b in xs if b and math.isfinite(a) and math.isfinite(b) and b > 0)
md.append(f"- held-out / predicted TTS, median over labels: primary {rat(pr):.3f}, secondary {rat(sr):.3f}")
# ---------------- power and energy ----------------
if os.path.exists('results/GSW/marks.txt'):
    marks = {}
    for l in open('results/GSW/marks.txt'):
        ts, name = l.split(); marks[name] = datetime.datetime.strptime(ts, '%Y-%m-%dT%H:%M:%S.%fZ')
    rows = []
    for l in open('results/GSW/power_samples.csv'):
        f = [x.strip() for x in l.split(',')]
        try: rows.append((datetime.datetime.strptime(f[0], '%Y/%m/%d %H:%M:%S.%f'), float(f[3].split()[0]), float(f[5].split()[0])))
        except Exception: pass
    off = rows[0][0] - marks['idle_start']
    def win(a, b, trim):
        w = [(p, m) for t, p, m in rows if marks[a] + off + datetime.timedelta(seconds=trim) <= t <= marks[b] + off - datetime.timedelta(seconds=trim)]
        return (sum(p for p, m in w) / len(w), sum(m for p, m in w) / len(w), len(w)) if w else (None, None, 0)
    P = {}
    for l in open('results/GSW/burns.tsv'):
        combo = l.split('\t')[0]; P[combo] = win(f'burn_{combo}_start', f'burn_{combo}_end', 5.0)
    idle = win('idle_start', min((k for k in marks if k.startswith('burn_') and k.endswith('_start')), key=lambda k: marks[k]), 2.0)
    out['power'] = dict(burns={k: dict(board_W=v[0], module_W=v[1], samples=v[2]) for k, v in P.items()}, idle_board_W=idle[0], idle_module_W=idle[1])
    md.append('\n### GPU board power at the new schedules (nvidia-smi on the UUID, 200 ms; 40 s per combination)\n\n| Combination (B, variant, class) | samples | board W | module W |\n|---|---:|---:|---:|')
    for k, v in P.items(): md.append(f"| {k} | {v[2]} | {fmt(v[0], 1)} | {fmt(v[1], 1)} |")
    md.append(f"| idle | {idle[2]} | {fmt(idle[0], 1)} | {fmt(idle[1], 1)} |")
    v80p = json.load(open('data/v80_a3/power_summary.json'))
    ph = v80p.get('phases', v80p)
    gl = [ph[k]['board_W_mean'] for k in ph if k.startswith('load_') and k != 'load_X5']
    v80W = sum(gl) / len(gl); out['v80_gset_board_W'] = v80W
    md.append(f"\nEnergy to solution, geometric means over the 51 instances (GPU: board power of the cohort's combination × TTS; V80: mean A3 board power of its G-set loads, {v80W:.1f} W, × TTS):\n")
    md.append('| Grid | Rule | GPU primary (mJ) | V80 primary (mJ) | GPU secondary (mJ) | V80 secondary A3 (mJ) |\n|---|---|---:|---:|---:|---:|')
    en = {}
    for grid in ('original', 'extended'):
        for r in RULES:
            def pw(k):
                e = est[G[k]['cohort']]; c = f"B{e['B']}_{e['variant']}_{'ons' if r == 'Onsager-online' else 'const'}"
                return P[c][0]
            gpe = gmean(pw((g, r, grid, 'primary')) * G[(g, r, grid, 'primary')]['primary'] for g in insts)
            gse = gmean(pw((g, r, grid, 'secondary')) * G[(g, r, grid, 'secondary')]['secondary'] for g in insts)
            vpe = gmean(v80W * V[(g, r, grid)]['primary'] for g in insts); vse = gmean(v80W * V[(g, r, grid)]['secondary'] for g in insts)
            en[f'{grid}|{r}'] = dict(gpu_primary_mJ=gpe, v80_primary_mJ=vpe, gpu_secondary_mJ=gse, v80_secondary_mJ=vse)
            md.append(f"| {grid} | {r} | {fmt(gpe, 1)} | {fmt(vpe, 2)} | {fmt(gse, 2)} | {fmt(vse, 3)} |")
    out['energy'] = en
json.dump(out, open('results/GS_tables.json', 'w'), indent=1, default=lambda x: None)
open('results/GS_tables.md', 'w').write('\n'.join(md) + '\n')
print('\n'.join(md))
