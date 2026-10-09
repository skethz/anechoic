#!/usr/bin/env python3
"""Collect the bias study's phases into results/tables.json and results/tables.md (run on gpu-host from the project root).
Estimators come from src/analyze.py (unchanged from the v2 / multibit studies)."""
import json, glob, math, os, sys, datetime, statistics
sys.path.insert(0, 'src')
import analyze as A

def fmt(x, d=4):
    if x is None: return '—'
    if isinstance(x, float) and math.isinf(x): return '∞'
    if isinstance(x, float) and math.isnan(x): return '—'
    return f'{x:.{d}f}'

def ci(c, d=4): return f'[{fmt(c[0], d)}, {fmt(c[1], d)}]'

def trials(prefix): return [json.loads(l) for l in open(prefix + '.trials.jsonl')]

def gmean(xs):
    xs = [x for x in xs if x is not None and x > 0 and not math.isinf(x)]
    return math.exp(sum(math.log(x) for x in xs) / len(xs)) if xs else None

out, md = {}, []
# ---------------- V ----------------
if os.path.exists('results/V/validation_summary.json'):
    vs = json.load(open('results/V/validation_summary.json'))
    out['validation'] = {k: vs[k] for k in vs if k not in ('rows', 'mb_rows')}
    groups = {}
    for r in vs['rows']:
        name = os.path.basename(r['prefix'])
        g = 'K2000 b=0' if name.startswith('k2000_') else 'K2000 zero-bias file' if name.startswith('k2000z_') else 'K2000 --wide 1' if name.startswith('k2000w_') else name.split('_')[0] + '_' + name.split('_')[1]
        if 'forcedwide' in name: g += ' --wide 1'
        x = groups.setdefault(g, dict(runs=0, trials=0, exact=0, N=r.get('N')))
        x['runs'] += 1; x['trials'] += r['trials']; x['exact'] += r['exact']
    out['validation']['groups'] = groups
    md.append(f"## Verification against run_trial_bias\n\n{vs['exact']}/{vs['trials']} trials bit-exact in {vs['runs']} runs (final spins, Σ s h, Σ b s, energy, flips, every n_lin); "
              f"{vs['mb_identical_trials']}/{vs['mb_identity_trials']} K2000 trials in {vs['mb_identity_runs']} runs identical to the multibit kernel's raw files; "
              f"device/host mismatches {vs['device_host_mismatches']}; final-field mismatches {vs['field_mismatches']}; WIDE runs {vs['wide_runs']}.\n")
    md.append('| Group | n | Runs | Trials | Exact |\n|---|---:|---:|---:|---:|')
    for g, x in groups.items(): md.append(f"| {g} | {x['N']} | {x['runs']} | {x['trials']} | {x['exact']} |")
# ---------------- T ----------------
if glob.glob('results/T/*.summary.json'):
    def t_of(pre): return json.load(open(pre + '.summary.json'))['device_ms_mean']
    T = {}
    md.append('\n## Step time\n\n### T1: bias binary (no bias) vs multibit binary, K2000, interleaved, 3 repetitions × 16 batches\n')
    md.append('| Config | S | B | bias t_batch (ms) | multibit t_batch (ms) | ratio (mean) | ratio range over reps |\n|---|---:|---:|---:|---:|---:|---|')
    S_of = dict(X5=280, X4=760, O4=360, O1=960)
    for c in ['X5', 'X4', 'O4', 'O1']:
        for B in [1, 60, 120, 240]:
            a = [t_of(f'results/T/t1_bias_{c}_b{B}_r{r}') for r in (1, 2, 3) if os.path.exists(f'results/T/t1_bias_{c}_b{B}_r{r}.summary.json')]
            b = [t_of(f'results/T/t1_mb_{c}_b{B}_r{r}') for r in (1, 2, 3) if os.path.exists(f'results/T/t1_mb_{c}_b{B}_r{r}.summary.json')]
            if not a or not b: continue
            rr = [x / y for x, y in zip(a, b)]
            T[f't1_{c}_b{B}'] = dict(bias=a, mb=b, ratio=statistics.mean(a) / statistics.mean(b), ratio_min=min(rr), ratio_max=max(rr))
            md.append(f"| {c} | {S_of[c]} | {B} | {fmt(statistics.mean(a))} | {fmt(statistics.mean(b))} | {fmt(statistics.mean(a) / statistics.mean(b), 4)} | {fmt(min(rr), 4)}–{fmt(max(rr), 4)} |")
    md.append('\nPer-step slopes (µs/step) from the two schedule lengths of a class:\n\n| Class | B | bias | multibit | ratio |\n|---|---:|---:|---:|---:|')
    for cls, (c1, c2) in [('const (X5→X4)', ('X5', 'X4')), ('Onsager (O4→O1)', ('O4', 'O1'))]:
        for B in [1, 60, 120, 240]:
            k1, k2 = f't1_{c1}_b{B}', f't1_{c2}_b{B}'
            if k1 not in T or k2 not in T: continue
            sb = 1000 * (statistics.mean(T[k2]['bias']) - statistics.mean(T[k1]['bias'])) / (S_of[c2] - S_of[c1])
            sm = 1000 * (statistics.mean(T[k2]['mb']) - statistics.mean(T[k1]['mb'])) / (S_of[c2] - S_of[c1])
            T[f'slope_{c1}{c2}_b{B}'] = dict(bias_us=sb, mb_us=sm, ratio=sb / sm)
            md.append(f"| {cls} | {B} | {fmt(sb, 4)} | {fmt(sm, 4)} | {fmt(sb / sm, 4)} |")
    md.append('\n### T2: with bias vs all-zero bias file (same binary, J, trial ids)\n\n| Instance | Mode | bias t_batch (ms) | zero-bias t_batch (ms) | ratio |\n|---|---|---:|---:|---:|')
    for inst in ['rb2_n1000', 'tsp40']:
        for mdn in ['on', 'pl']:
            a = [t_of(f'results/T/t2_bias_{inst}_{mdn}_r{r}') for r in (1, 2, 3) if os.path.exists(f'results/T/t2_bias_{inst}_{mdn}_r{r}.summary.json')]
            b = [t_of(f'results/T/t2_zero_{inst}_{mdn}_r{r}') for r in (1, 2, 3) if os.path.exists(f'results/T/t2_zero_{inst}_{mdn}_r{r}.summary.json')]
            if not a or not b: continue
            T[f't2_{inst}_{mdn}'] = dict(bias=a, zero=b, ratio=statistics.mean(a) / statistics.mean(b))
            md.append(f"| {inst} | {mdn} | {fmt(statistics.mean(a))} | {fmt(statistics.mean(b))} | {fmt(statistics.mean(a) / statistics.mean(b), 4)} |")
    md.append('\n### T3: clamped (WIDE) vs plain kernel, K2000\n\n| Config | B | WIDE t_batch (ms) | plain t_batch (ms) | ratio |\n|---|---:|---:|---:|---:|')
    for c in ['X5', 'O1']:
        for B in [1, 120, 240]:
            a = [t_of(f'results/T/t3_wide_{c}_b{B}_r{r}') for r in (1, 2, 3) if os.path.exists(f'results/T/t3_wide_{c}_b{B}_r{r}.summary.json')]
            b = [t_of(f'results/T/t3_plain_{c}_b{B}_r{r}') for r in (1, 2, 3) if os.path.exists(f'results/T/t3_plain_{c}_b{B}_r{r}.summary.json')]
            if not a or not b: continue
            T[f't3_{c}_b{B}'] = dict(wide=a, plain=b, ratio=statistics.mean(a) / statistics.mean(b))
            md.append(f"| {c} | {B} | {fmt(statistics.mean(a))} | {fmt(statistics.mean(b))} | {fmt(statistics.mean(a) / statistics.mean(b), 4)} |")
    out['T'] = T
    if os.path.exists('results/T/sass_compare.json'):
        sc = json.load(open('results/T/sass_compare.json')); out['T_sass'] = sc
        md.append('\n### Static SASS comparison (bias kernel vs multibit kernel)\n\n| Kernel | multibit instructions | bias instructions | opcode diff (inserted / deleted) |\n|---|---:|---:|---|')
        for k, v in sc.items(): md.append(f"| {k} | {v['mb']} | {v['bias']} | {v['inserted']} / {v['deleted']} |")
# ---------------- BV / BT (Amendment 1) ----------------
if os.path.exists('results/BV/summary.json'):
    bv = json.load(open('results/BV/summary.json')); out['BV'] = bv
    md.append(f"\n## Binary sca_gpu_bias_tg (Amendment 1)\n\nVerification: {bv['exact']}/{bv['trials']} trials bit-exact in {bv['runs']} runs ({bv['tables_global_runs']} with global tables, "
              f"including S = 5000-8192); {bv['mb_identical_trials']}/{bv['mb_identity_trials']} K2000 trials identical to the multibit kernel; device/host mismatches {bv['device_host_mismatches']}; field mismatches {bv['field_mismatches']}.\n")
if glob.glob('results/BT/*.summary.json'):
    def tt(pre): return json.load(open(pre + '.summary.json'))['device_ms_mean']
    BT = {}
    md.append('Global vs shared-memory step tables (3 interleaved repetitions × 16 batches):\n\n| Run | B | global t_batch (ms) | smem t_batch (ms) | ratio |\n|---|---:|---:|---:|---:|')
    for name in ['k2000_X5', 'k2000_O1', 'tsp40_on', 'tsp40_pl']:
        for B in [1, 60, 120, 240]:
            a = [tt(f'results/BT/{name}_b{B}_global_r{r}') for r in (1, 2, 3) if os.path.exists(f'results/BT/{name}_b{B}_global_r{r}.summary.json')]
            b = [tt(f'results/BT/{name}_b{B}_smem_r{r}') for r in (1, 2, 3) if os.path.exists(f'results/BT/{name}_b{B}_smem_r{r}.summary.json')]
            if not a or not b: continue
            BT[f'{name}_b{B}'] = dict(glob=a, smem=b, ratio=statistics.mean(a) / statistics.mean(b))
            md.append(f"| {name} | {B} | {fmt(statistics.mean(a))} | {fmt(statistics.mean(b))} | {fmt(statistics.mean(a) / statistics.mean(b), 4)} |")
    out['BT'] = BT
# ---------------- GV / G ----------------
CLASS = {}
for g in list(range(1, 11)) + list(range(22, 32)) + list(range(43, 48)): CLASS[f'G{g}'] = 'random'
for g in list(range(11, 14)) + list(range(32, 35)): CLASS[f'G{g}'] = 'toroidal'
for g in list(range(14, 22)) + list(range(35, 43)) + list(range(51, 55)): CLASS[f'G{g}'] = 'planar'
if os.path.exists('results/GV/summary.json'):
    out['G_validation'] = json.load(open('results/GV/summary.json'))
    gv = out['G_validation']
    md.append(f"\n## G-set\n\nValidation (before counting): {gv['exact']}/{gv['trials']} trials bit-exact against run_trial_bias(n = N, b = 0) in {gv['runs']} runs; device/host mismatches {gv['device_host_mismatches']}.\n")
if os.path.exists('results/G/plan.json') and glob.glob('results/G/*_b*.summary.json'):
    plan = json.load(open('results/G/plan.json'))
    weights = {}
    for g in set(c['instance'] for c in plan):
        neg = False
        with open(f'data/gset/{g}') as f:
            f.readline()
            for l in f:
                if int(l.split()[2]) < 0: neg = True; break
        weights[g] = '±1' if neg else '+1'
    G = {}
    for c in plan:
        for B in ([60, 120, 240] if c['kind'] == 'E12' else [1]):
            pre = f"results/G/{c['tag']}_b{B}"
            if not os.path.exists(pre + '.summary.json'): continue
            e = A.estimate(pre); s = json.load(open(pre + '.summary.json')); tr = trials(pre)
            cuts = [x['cut'] for x in tr]
            e.update(instance=c['instance'], N=c['N'], target=c['target'], BKV=c['BKV'], kind=c['kind'], rule=c['rule'], cls=c['cls'],
                     graph_class=CLASS[c['instance']], weights=weights[c['instance']], variant=f"R{s['r']}G{s['groups']}", mean_cut=sum(cuts) / len(cuts),
                     best_cut=max(cuts), study_tts_ms=c['study_tts_ms'], study_p=c['study_p'], mismatches=s['mismatches'])
            G[e['prefix']] = e
    out['G'] = G
    rules = ['SCA', 'TEC', 'Onsager-kT', 'Onsager-online']
    # per instance: E12 primary TTS at each B per rule (best B), E1 at B = 1
    md.append('### Per instance: primary TTS99 (ms) at the study\'s E12 configuration (best of B = 60/120/240, B in parentheses) and single-chain TTS99 at the E1 configuration\n')
    md.append('| Instance | N | class | w | target | ' + ' | '.join(f'{r} E12' for r in rules) + ' | ' + ' | '.join(f'{r} E1 B=1' for r in rules) + ' | best cut (all runs) |')
    md.append('|---|---:|---|---|---:|' + '---:|' * (2 * len(rules)) + '---:|')
    insts = sorted(set(c['instance'] for c in plan), key=lambda g: int(g[1:]))
    best_per_inst = {}
    for g in insts:
        row = []; allbest = -10**9; bestp = None
        for r in rules:
            tag = f"{g}_{r.replace('-', '')}_E12"
            cand = [(G[f'{tag}_b{B}']['tts_primary_ms'], B) for B in (60, 120, 240) if f'{tag}_b{B}' in G]
            for B in (60, 120, 240):
                if f'{tag}_b{B}' in G: allbest = max(allbest, G[f'{tag}_b{B}']['best_cut'])
            if cand:
                t, B = min(cand); row.append(f"{fmt(t, 3)} ({B})")
                if bestp is None or t < bestp[0]: bestp = (t, r, B)
            else: row.append('—')
        for r in rules:
            k = f"{g}_{r.replace('-', '')}_E1_b1"
            if k in G: row.append(fmt(G[k]['tts_primary_ms'], 3)); allbest = max(allbest, G[k]['best_cut'])
            else: row.append('—')
        c0 = next(c for c in plan if c['instance'] == g)
        best_per_inst[g] = dict(best_primary=bestp, best_cut=allbest, N=c0['N'], target=c0['target'], BKV=c0['BKV'], graph_class=CLASS[g], weights=weights[g])
        md.append(f"| {g} | {c0['N']} | {CLASS[g]} | {weights[g]} | {c0['target']} | " + ' | '.join(row) + f" | {allbest} (BKV {c0['BKV']}) |")
    out['G_best_per_instance'] = best_per_inst
    # aggregates by class
    md.append('\n### By graph class: geometric mean over instances of the primary TTS99 (ms), E12 configuration at B = 120, and fraction of cells reaching the target at least once\n')
    md.append('| Class | N | w | instances | ' + ' | '.join(rules) + ' | best rule per instance (geomean) |\n|---|---:|---|---:|' + '---:|' * (len(rules) + 1))
    keys = sorted(set((best_per_inst[g]['graph_class'], best_per_inst[g]['N'], best_per_inst[g]['weights']) for g in insts), key=lambda k: (k[1], k[0], k[2]))
    agg = {}
    for (cl, N, w) in keys:
        gs = [g for g in insts if (best_per_inst[g]['graph_class'], best_per_inst[g]['N'], best_per_inst[g]['weights']) == (cl, N, w)]
        cells = []
        for r in rules:
            vals = [G.get(f"{g}_{r.replace('-', '')}_E12_b120", {}).get('tts_primary_ms') for g in gs]
            reach = sum(1 for g in gs if G.get(f"{g}_{r.replace('-', '')}_E12_b120", {}).get('successes', 0) > 0)
            gm = gmean(vals) if reach == len(gs) else None
            cells.append(f"{fmt(gm, 3)} ({reach}/{len(gs)})")
        bp = [best_per_inst[g]['best_primary'][0] if best_per_inst[g]['best_primary'] else None for g in gs]
        agg[f'{cl}_{N}_{w}'] = dict(instances=gs, best_primary_geomean=gmean(bp))
        md.append(f"| {cl} | {N} | {w} | {len(gs)} | " + ' | '.join(cells) + f" | {fmt(gmean(bp), 3)} |")
    out['G_classes'] = agg
    nm = sum(e['mismatches'] for e in G.values())
    md.append(f"\nAll {len(G)} counted G-set runs: device/host cut mismatches {nm}.")
    # transfer check: GPU p (pooled over B for E12) vs the study's float-model p for the same configuration
    import math as _m
    rows = []
    for c in plan:
        Bs = [60, 120, 240] if c['kind'] == 'E12' else [1]
        ks = [G[f"{c['tag']}_b{B}"] for B in Bs if f"{c['tag']}_b{B}" in G]
        if not ks: continue
        k = sum(e['successes'] for e in ks); n = sum(e['trials'] for e in ks)
        # two-proportion z test; the study's cohort is taken as 1,024 trials (p * 1024 is an integer in every cell)
        n2 = 1024; k2 = round(c['study_p'] * n2); pp = (k + k2) / (n + n2)
        se = _m.sqrt(pp * (1 - pp) * (1 / n + 1 / n2)) if 0 < pp < 1 else 0.0
        z = (k / n - k2 / n2) / se if se > 0 else 0.0
        rows.append(dict(tag=c['tag'], kind=c['kind'], p_gpu=k / n, p_study=c['study_p'], n=n, z=z, inside=abs(z) < 1.959963984540054, diff=k / n - c['study_p']))
    out['G_transfer'] = rows
    for kind in ['E12', 'E1']:
        r = [x for x in rows if x['kind'] == kind]
        if not r: continue
        d = sorted(abs(x['diff']) for x in r)
        md.append(f"\nSuccess probability, GPU vs the study's float model at the same {kind} configuration ({len(r)} cells): median |Δp| = {d[len(d)//2]:.4f}, "
                  f"max {d[-1]:.4f}, mean Δp = {sum(x['diff'] for x in r)/len(r):+.4f}; two-proportion test (study cohort 1,024 trials): |z| < 1.96 in {sum(x['inside'] for x in r)}/{len(r)} cells.")
    # per instance: best secondary and the study's modelled V80 times (12-engine and single-engine), best over rules
    md.append('\n### Per instance: best over rules, GPU against the study\'s modelled V80 times (study numbers from selected_configs.json)\n')
    md.append('| Instance | class | N | GPU best primary (ms) [rule, B] | GPU best secondary (ms) [rule, B] | study V80 E12 model best (ms) [rule] | GPU best B=1 (ms) [rule] | study V80 E1 model best (ms) [rule] | mean cut / BKV (best rule E12, B=240) |')
    md.append('|---|---|---:|---|---|---|---|---|---:|')
    cmp_rows = {}
    for g in insts:
        cand_p, cand_s, cand_1, st12, st1, mq = [], [], [], [], [], []
        for r in rules:
            tag = f"{g}_{r.replace('-', '')}"
            for B in (60, 120, 240):
                e = G.get(f'{tag}_E12_b{B}')
                if e: cand_p.append((e['tts_primary_ms'], r, B)); cand_s.append((e['tts_secondary_ms'], r, B))
            e1 = G.get(f'{tag}_E1_b1')
            if e1: cand_1.append((e1['tts_primary_ms'], r))
            c12 = next(c for c in plan if c['tag'] == f'{tag}_E12'); c1 = next(c for c in plan if c['tag'] == f'{tag}_E1')
            st12.append((c12['study_tts_ms'], r)); st1.append((c1['study_tts_ms'], r))
            e240 = G.get(f'{tag}_E12_b240')
            if e240: mq.append((e240['tts_primary_ms'], e240['mean_cut'] / e240['BKV'], r))
        bp, bs, b1, s12, s1 = min(cand_p), min(cand_s), min(cand_1), min(st12), min(st1)
        q = min(mq)[1] if mq else None
        cmp_rows[g] = dict(gpu_primary=bp, gpu_secondary=bs, gpu_b1=b1, study_e12=s12, study_e1=s1, mean_cut_over_bkv=q)
        md.append(f"| {g} | {CLASS[g]} {weights[g]} | {next(c['N'] for c in plan if c['instance'] == g)} | {fmt(bp[0], 3)} [{bp[1]}, {bp[2]}] | {fmt(bs[0], 4)} [{bs[1]}, {bs[2]}] | "
                  f"{fmt(s12[0], 4)} [{s12[1]}] | {fmt(b1[0], 3)} [{b1[1]}] | {fmt(s1[0], 4)} [{s1[1]}] | {fmt(q, 4)} |")
    out['G_compare'] = cmp_rows
    md.append('\n### By class: geometric means over instances of the per-instance bests (ms)\n\n| Class | N | w | n | GPU primary | GPU secondary | study V80 E12 model | GPU B=1 | study V80 E1 model |\n|---|---:|---|---:|---:|---:|---:|---:|---:|')
    for (cl, N, w) in keys:
        gs = [g for g in insts if (best_per_inst[g]['graph_class'], best_per_inst[g]['N'], best_per_inst[g]['weights']) == (cl, N, w)]
        f = lambda key: gmean([cmp_rows[g][key][0] for g in gs])
        md.append(f"| {cl} | {N} | {w} | {len(gs)} | {fmt(f('gpu_primary'), 3)} | {fmt(f('gpu_secondary'), 4)} | {fmt(f('study_e12'), 4)} | {fmt(f('gpu_b1'), 3)} | {fmt(f('study_e1'), 4)} |")
    f = lambda key: gmean([cmp_rows[g][key][0] for g in insts])
    md.append(f"| all | | | {len(insts)} | {fmt(f('gpu_primary'), 3)} | {fmt(f('gpu_secondary'), 4)} | {fmt(f('study_e12'), 4)} | {fmt(f('gpu_b1'), 3)} | {fmt(f('study_e1'), 4)} |")
    # energy to solution from the measured board power of the G-set burn (phase W: G22 Onsager-online E12, B = 120)
    if os.path.exists('results/W/marks.txt'):
        Pw = None
        try:
            import csv
            mk = {l.split()[1]: datetime.datetime.strptime(l.split()[0], '%Y-%m-%dT%H:%M:%S.%fZ') for l in open('results/W/marks.txt')}
            rows_ = []
            for l in open('results/W/power_samples.csv'):
                f = [x.strip() for x in l.split(',')]
                try: rows_.append((datetime.datetime.strptime(f[0], '%Y/%m/%d %H:%M:%S.%f'), float(f[3].split()[0])))
                except Exception: pass
            off = rows_[0][0] - mk['idle_start']
            win = [v for t, v in rows_ if mk['burn_g_start'] + off + datetime.timedelta(seconds=5) <= t <= mk['burn_g_end'] + off - datetime.timedelta(seconds=5)]
            idl = [v for t, v in rows_ if mk['idle_start'] + off + datetime.timedelta(seconds=2) <= t <= mk['burn_g_start'] + off - datetime.timedelta(seconds=2)]
            Pw, Pi = sum(win) / len(win), sum(idl) / len(idl)
        except Exception:
            Pw = None
        if Pw:
            fP = lambda key: gmean([cmp_rows[g][key][0] for g in insts])
            out['G_energy'] = dict(board_W=Pw, idle_W=Pi, primary_mJ=Pw * fP('gpu_primary'), primary_idle_subtracted_mJ=(Pw - Pi) * fP('gpu_primary'),
                                   secondary_mJ=Pw * fP('gpu_secondary'), single_chain_mJ=Pw * fP('gpu_b1'))
            md.append(f"\nEnergy to solution (board power {Pw:.1f} W measured on the G22 burn, idle {Pi:.1f} W; geometric means over the 51 instances of the per-instance best): "
                      f"primary {Pw * fP('gpu_primary'):.0f} mJ ({(Pw - Pi) * fP('gpu_primary'):.0f} mJ idle-subtracted), secondary {Pw * fP('gpu_secondary'):.1f} mJ, "
                      f"single chain {Pw * fP('gpu_b1'):.0f} mJ.")
    nb = sum(1 for g in insts if best_per_inst[g]['best_cut'] >= best_per_inst[g]['BKV'])
    md.append(f"\nBest cut over all GPU runs equals the best-known value on {nb}/{len(insts)} instances (it never exceeds it).")
# ---------------- B (Amendment 3): graph partitioning and TSP ----------------
if os.path.exists('data/hw_export/import_check.json') and glob.glob('results/B/*.eval.json'):
    ic = json.load(open('data/hw_export/import_check.json')); out['B_import'] = ic
    plan_b = [l.rstrip('\n').split('\t') for l in open('data/hw_export/plan_b.tsv')]
    if os.path.exists('results/BV2/summary.json'):
        bv2 = json.load(open('results/BV2/summary.json')); out['BV2'] = bv2
        md.append(f"\n## Graph partitioning and TSP (K = 4, ReAIM study export)\n\nImport: {len(ic['files'])} files hash-checked; {len(ic['instances'])} instances with full mapping checks "
                  f"(matrix and bias rebuilt from the original problem: {sum(v['mapping_check_J'] and v['mapping_check_b'] for v in ic['instances'].values())}/{len(ic['instances'])}); "
                  f"{sum(x['tables_identical'] for x in ic['schedules'])}/{len(ic['schedules'])} schedules with tables identical to mb_common.hpp::tables. "
                  f"Validation: {bv2['exact']}/{bv2['trials']} traced trials bit-exact against run_trial_bias in {bv2['runs']} runs.\n")
    B = {}
    md.append('| Instance | Rule | S | global tables | replication of the study (1024 trials): feasible / P(target) / mean quality, GPU = model? | B=60: feasibility | P(target) [95%] | mean quality | best value (ref) | TTS99 primary (ms) [95%] | TTS99 secondary (ms) | B=1 TTS99 (ms) |')
    md.append('|---|---|---:|---|---|---:|---|---:|---|---|---:|---:|')
    for tag, inst, kind, N, rule, S, seed, jf, bf, spec, fl in plan_b:
        e60 = json.load(open(f'results/B/{tag}_b60.eval.json')) if os.path.exists(f'results/B/{tag}_b60.eval.json') else None
        e1 = json.load(open(f'results/B/{tag}_b1.eval.json')) if os.path.exists(f'results/B/{tag}_b1.eval.json') else None
        er = json.load(open(f'results/BR/{tag}.eval.json')) if os.path.exists(f'results/BR/{tag}.eval.json') else None
        sch = next(x for x in ic['schedules'] if x['instance'] == inst and x['rule'] == rule); ps = sch.get('precision_study') or {}
        rep = None
        if er and ps:
            rep = dict(gpu=dict(p_feasible=er['feasibility'], p_target=er['p'], mean_quality=er['mean_quality']),
                       model=dict(p_feasible=ps.get('p_feasible'), p_target=ps.get('p_target'), mean_quality=ps.get('mean_quality')))
            rep['identical'] = (abs(er['feasibility'] - ps.get('p_feasible', -1)) < 1e-12 and abs(er['p'] - ps.get('p_target', -1)) < 1e-12
                                and abs(er['mean_quality'] - ps.get('mean_quality', -1)) < 1e-9)
        B[tag] = dict(instance=inst, kind=kind, rule=rule, S=int(S), b60=e60, b1=e1, replication=rep, tables_global=sch['tables_global'])
        ref = json.load(open(spec))['ref']
        reptxt = '—' if not rep else f"{rep['gpu']['p_feasible']:.4f} / {rep['gpu']['p_target']:.4f} / {rep['gpu']['mean_quality']:.4f}, {'yes' if rep['identical'] else 'NO'}"
        if e60:
            md.append(f"| {inst} | {rule} | {S} | {'yes' if sch['tables_global'] else 'no'} | {reptxt} | {fmt(e60['feasibility'], 3)} | {fmt(e60['p'], 4)} {ci(e60['p_wilson95'], 4)} | "
                      f"{fmt(e60['mean_quality'], 4)} | {e60['best_value']} ({ref}) | {fmt(e60['tts_primary_ms'], 3)} {ci(e60['tts_primary_ci_ms'], 3)} | "
                      f"{fmt(e60['tts_secondary_ms'], 4)} | {fmt(e1['tts_primary_ms'], 2) if e1 else '—'} |")
    out['B'] = B
    nrep = [v['replication']['identical'] for v in B.values() if v['replication']]
    md.append(f"\nReplication of the study's bit-exact model on its own trials: {sum(nrep)}/{len(nrep)} cells identical (feasibility, P(target), mean quality).")
    diffc = [k for k, v in B.items() if v['replication'] and not v['replication']['identical']]
    if diffc:
        md.append('Cells that differ: ' + ', '.join(f"{k} (GPU feasibility {B[k]['replication']['gpu']['p_feasible']:.3f} vs exported {B[k]['replication']['model']['p_feasible']:.3f})" for k in diffc) + '.')
    if os.path.exists('results/BR_diag/check.jsonl'):
        dg = [json.loads(l) for l in open('results/BR_diag/check.jsonl')]
        md.append(f"Post-hoc diagnostic (not pre-registered): on trials 0-31 of each differing cell the GPU is bit-exact with run_trial_bias in "
                  f"{sum(x['exact'] for x in dg)}/{sum(x['trials'] for x in dg)} trials ({len(dg)} cells), so the golden reference itself gives the GPU's values there.")
    # per problem: best rule per instance by B=60 mean quality
    md.append('\n### Best rule per instance (B = 60, by mean quality)\n\n| Instance | kind | rule | S | feasibility | P(target) | mean quality | best value (ref) | TTS99 primary (ms) | TTS99 secondary (ms) |\n|---|---|---|---:|---:|---:|---:|---|---:|---:|')
    for inst in dict.fromkeys(v['instance'] for v in B.values()):
        cs = [(k, v) for k, v in B.items() if v['instance'] == inst and v['b60']]
        if not cs: continue
        k, v = max(cs, key=lambda kv: kv[1]['b60']['mean_quality']); e = v['b60']; ref = json.load(open(f"data/hw_export/{inst}_K4.spec.json"))['ref']
        md.append(f"| {inst} | {v['kind']} | {v['rule']} | {v['S']} | {fmt(e['feasibility'], 3)} | {fmt(e['p'], 4)} | {fmt(e['mean_quality'], 4)} | {e['best_value']} ({ref}) | {fmt(e['tts_primary_ms'], 3)} | {fmt(e['tts_secondary_ms'], 4)} |")
# ---------------- W ----------------
def power_table(WD, title, labels):
    marks = {}
    for l in open(f'results/{WD}/marks.txt'):
        ts, name = l.split(); marks[name] = datetime.datetime.strptime(ts, '%Y-%m-%dT%H:%M:%S.%fZ')
    rows = []
    for l in open(f'results/{WD}/power_samples.csv'):
        f = [x.strip() for x in l.split(',')]
        try:
            ts = datetime.datetime.strptime(f[0], '%Y/%m/%d %H:%M:%S.%f')
            rows.append((ts, [float(x.split()[0]) if x.split()[0] not in ('[N/A]', 'N/A') else math.nan for x in f[1:8]]))
        except Exception:
            pass
    off = rows[0][0] - marks['idle_start']
    def window(a, b, trim=2.0):
        return [v for ts, v in rows if (marks[a] + off + datetime.timedelta(seconds=trim)) <= ts <= (marks[b] + off - datetime.timedelta(seconds=trim))]
    names = ['util', 'mem_mib', 'power_draw_w', 'power_draw_instant_w', 'module_power_avg_w', 'module_power_instant_w', 'sm_clock_mhz']
    def stats(ws): return {n: (sum(v[i] for v in ws) / len(ws) if ws else math.nan) for i, n in enumerate(names)} | {'samples': len(ws)}
    order = [k for k in ['idle_start', 'burn_g_start', 'burn_g_end', 'burn_b_start', 'burn_b_end', 'idle_end'] if k in marks]
    W = {}
    md.append(f'\n## {title}\n\n| Window | Samples | util % | board W (power.draw) | module W (avg) | SM MHz |\n|---|---:|---:|---:|---:|---:|')
    for a, b in zip(order[:-1], order[1:]):
        st = stats(window(a, b, 5.0 if a.startswith('burn') else 2.0)); W[f'{a}->{b}'] = st
        md.append(f"| {labels.get(a, a)} | {st['samples']} | {fmt(st['util'], 1)} | {fmt(st['power_draw_w'], 1)} | {fmt(st['module_power_avg_w'], 1)} | {fmt(st['sm_clock_mhz'], 0)} |")
    for bj in sorted(glob.glob(f'results/{WD}/*.burn.json')):
        d = json.load(open(bj)); W[os.path.basename(bj)] = d
        md.append(f"\n{os.path.basename(bj)}: {d['batches']} batches of B = {d['chains']} in {d['burn_seconds_wall']:.1f} s, {d['device_ms_per_batch']:.4f} ms per batch.")
    out[WD] = W
if os.path.exists('results/W/power_samples.csv'):
    power_table('W', 'Power, G-set and biased TSP40 (GPU1 by UUID, 200 ms sampling)', {'idle_start': 'idle', 'burn_g_start': 'G22 Onsager-online E12, B = 120', 'burn_g_end': 'idle', 'burn_b_start': 'tsp40 K = 4 biased, B = 60', 'burn_b_end': 'idle'})
if os.path.exists('results/BW/power_samples.csv'):
    power_table('BW', 'Power, GP/TSP (GPU1 by UUID, 200 ms sampling)', {'idle_start': 'idle', 'burn_g_start': 'TSP burn', 'burn_g_end': 'idle', 'burn_b_start': 'GPP burn', 'burn_b_end': 'idle'})
json.dump(out, open('results/tables.json', 'w'), indent=1, default=float)
open('results/tables.md', 'w').write('\n'.join(md) + '\n')
print('\n'.join(md))
