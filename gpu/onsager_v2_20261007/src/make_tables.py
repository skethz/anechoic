#!/usr/bin/env python3
"""Collect all phases into results/tables.json and results/tables.md (run from the project root on gpu-host)."""
import json, glob, math, os, sys, datetime
sys.path.insert(0, 'src')
import analyze as A

V80 = {  # fpga/v80_sca/RESULTS_HW.md, 12 x v6.2 at 275 MHz (board_v62_e12_275mhz)
    'X5': dict(primary=0.0686, secondary=0.0638, p=0.338, t_round=0.0686),
    'X1': dict(primary=0.0771, secondary=0.0529, p=0.428, t_round=0.0771),
    'X4': dict(primary=0.1984, secondary=0.0254, p=0.950, t_round=0.1984),
    'O1': dict(primary=0.2352, secondary=0.0333, p=0.933, t_round=0.2352),
    'O4': dict(primary=0.0695, secondary=0.0756, p=0.326, t_round=0.0775),
}
ORDER = ['X5', 'X1', 'X4', 'X2', 'O1', 'O4', 'P4']


def mean_flips(prefix):
    n = 0; s = 0
    for l in open(prefix + '.trials.jsonl'):
        d = json.loads(l); s += d['flips']; n += 1
    return s / n


def fmt(x, d=4):
    if x is None: return '—'
    if isinstance(x, float) and math.isinf(x): return '∞'
    return f'{x:.{d}f}'


def ci(c, d=4):
    return f'[{fmt(c[0], d)}, {fmt(c[1], d)}]'


out = {}
md = []
# ---- V ----
vs = json.load(open('results/V/validation_summary.json'))
out['validation'] = vs
md.append(f"## Validation\n\n{vs['exact']}/{vs['trials']} trials bit-exact in {vs['runs']} runs (spins, cut, flips, every n_lin).\n")
md.append('| Run | Trials | Exact |\n|---|---:|---:|')
for r in vs['rows']:
    md.append(f"| {os.path.basename(r['prefix'])} | {r['trials']} | {r['exact']} |")
# ---- K ----
md.append('\n## Kernel variant selection (X5, 16 batches; 4 for L2)\n\n| Run | B | CS | G | t_batch (ms) | µs/step |\n|---|---:|---:|---:|---:|---:|')
out['K'] = {}
for f in sorted(glob.glob('results/K/*.summary.json')):
    e = A.estimate(f[:-13]); out['K'][e['prefix']] = e
    md.append(f"| {e['prefix']} | {e['B']} | {e['cs']} | {e['groups']} | {fmt(e['t_batch_ms'])} | {fmt(e['us_per_step'], 3)} |")
md.append('\nSelected: ' + json.dumps(json.load(open('results/K/selected_variants.json'))))
# ---- M ----
out['M'] = {}
md.append('\n## V80 configurations\n\n| Config | B | Batches | p [Wilson 95%] | P_batch | t_batch (ms) | µs/step | TTS99 primary (ms) [95%] | TTS99 secondary (ms) [95%] | V80 primary / secondary (ms) |')
md.append('|---|---:|---:|---|---:|---:|---:|---|---|---|')
for c in ORDER:
    for B in [1, 120, 240, 480]:
        pre = f'results/M/{c}_b{B}'
        if not os.path.exists(pre + '.summary.json'): continue
        e = A.estimate(pre); e['mean_flips'] = mean_flips(pre); out['M'][f'{c}_b{B}'] = e
        v = V80.get(c)
        vtxt = f"{v['primary']} / {v['secondary']}" if v and B in (120, 240) else ''
        md.append(f"| {c} | {B} | {e['batches']} | {fmt(e['p'], 3)} {ci(e['p_wilson95'], 3)} | {fmt(e['P_batch'], 3)} | {fmt(e['t_batch_ms'])} | {fmt(e['us_per_step'], 3)} | "
                  f"{fmt(e['tts_primary_ms'])} {ci(e['tts_primary_ci_ms'])} | {fmt(e['tts_secondary_ms'])} {ci(e['tts_secondary_ci_ms'])} | {vtxt} |")
# ---- P ----
if os.path.exists('results/P/pilot_cells.json'):
    cells = json.load(open('results/P/pilot_cells.json'))
    for c in cells: c['mean_flips'] = mean_flips('results/P/' + c['prefix'])
    out['P'] = cells
    # per-step fits per variant/mode: t = a + b S (+ c flips_per_trial)
    import numpy as np
    fits = {}
    for B in [120, 240]:
        for ons in [False, True]:
            sel = [c for c in cells if c['B'] == B and ((c['lam'] > 0) == ons)]
            if len(sel) < 4: continue
            S = np.array([c['steps'] for c in sel], float); t = np.array([c['t_batch_ms'] for c in sel]); F = np.array([c['mean_flips'] for c in sel])
            X1 = np.vstack([np.ones_like(S), S]).T; b1, res1, *_ = np.linalg.lstsq(X1, t, rcond=None)
            X2 = np.vstack([np.ones_like(S), S, F]).T; b2, *_ = np.linalg.lstsq(X2, t, rcond=None)
            r1 = t - X1 @ b1
            fits[f'B{B}_{"onsager" if ons else "const_corr"}'] = dict(n=len(sel), a_ms=b1[0], us_per_step=1000 * b1[1], max_abs_resid_ms=float(np.max(np.abs(r1))),
                                                                        with_flips=dict(a_ms=b2[0], us_per_step=1000 * b2[1], ns_per_flip_per_trial=1e6 * b2[2]))
    out['P_fits'] = fits
    md.append('\n## Per-step time fits over pilot cells (t_batch = a + b·S [+ c·mean flips per trial])\n\n| Cells | n | a (ms) | b (µs/step) | max |resid| (ms) | with flips: b (µs/step), c (ns per flip per trial) |\n|---|---:|---:|---:|---:|---|')
    for k, f in fits.items():
        md.append(f"| {k} | {f['n']} | {fmt(f['a_ms'])} | {fmt(f['us_per_step'], 4)} | {fmt(f['max_abs_resid_ms'])} | {fmt(f['with_flips']['us_per_step'], 4)}, {fmt(f['with_flips']['ns_per_flip_per_trial'], 3)} |")
    for w in ['primary', 'secondary']:
        top = sorted(cells, key=lambda c: c['score_' + w])[:8]
        md.append(f'\n### Pilot (exploratory) — best 8 cells by the pre-registered {w} score\n\n| Cell | Schedule | B | p | P_batch | t_batch (ms) | score (ms) | point TTS (ms) |\n|---|---|---:|---:|---:|---:|---:|---:|')
        for c in top:
            mode = f"TEC-T κ{c['kappa']:g}" if c['kappa'] > 0 else f"Onsager λ{c['lam']:g}"
            md.append(f"| {c['prefix']} | {mode}+ramp q{c['q']:g} T{c['t0']:g} S{c['steps']} | {c['B']} | {fmt(c['p'], 3)} | {fmt(c['P_batch'], 3)} | {fmt(c['t_batch_ms'])} | {fmt(c['score_' + w])} | {fmt(c['tts_' + w + '_ms'])} |")
# ---- H ----
if glob.glob('results/H/*.summary.json'):
    out['H'] = {}
    md.append('\n## Held-out (GPU-selected schedules)\n\n' + open('results/H/selected.txt').read().strip() + '\n')
    md.append('| Run | B | Batches | p [95%] | P_batch [95%] | t_batch (ms) | TTS99 primary (ms) [95%] | TTS99 secondary (ms) [95%] |\n|---|---:|---:|---|---|---:|---|---|')
    for f in sorted(glob.glob('results/H/*.summary.json')):
        e = A.estimate(f[:-13]); e['mean_flips'] = mean_flips(f[:-13]); out['H'][e['prefix']] = e
        md.append(f"| {e['prefix']} | {e['B']} | {e['batches']} | {fmt(e['p'], 3)} {ci(e['p_wilson95'], 3)} | {fmt(e['P_batch'], 3)} {ci(e['P_batch_wilson95'], 3)} | {fmt(e['t_batch_ms'])} | "
                  f"{fmt(e['tts_primary_ms'])} {ci(e['tts_primary_ci_ms'])} | {fmt(e['tts_secondary_ms'])} {ci(e['tts_secondary_ci_ms'])} |")
# ---- W ----
if os.path.exists('results/W/power_samples.csv'):
    marks = {}
    for l in open('results/W/marks.txt'):
        ts, name = l.split(); marks[name] = datetime.datetime.strptime(ts, '%Y-%m-%dT%H:%M:%S.%fZ')
    rows = []
    for l in open('results/W/power_samples.csv'):
        f = [x.strip() for x in l.split(',')]
        try:
            ts = datetime.datetime.strptime(f[0], '%Y/%m/%d %H:%M:%S.%f')
            vals = [float(x.split()[0]) if x.split()[0] not in ('[N/A]', 'N/A') else math.nan for x in f[1:8]]
            rows.append((ts, vals))
        except Exception:
            pass
    # nvidia-smi timestamps are local time; align by the offset between the first sample and the idle_start mark
    off = rows[0][0] - marks['idle_start']
    def window(a, b, trim=2.0):
        return [v for ts, v in rows if (marks[a] + off + datetime.timedelta(seconds=trim)) <= ts <= (marks[b] + off - datetime.timedelta(seconds=trim))]
    names = ['util', 'mem_mib', 'power_draw_w', 'power_draw_instant_w', 'module_power_avg_w', 'module_power_instant_w', 'sm_clock_mhz']
    def stats(ws):
        return {n: (sum(v[i] for v in ws) / len(ws) if ws else math.nan) for i, n in enumerate(names)} | {'samples': len(ws)}
    idle = stats(window('idle_start', 'burn_start')); burn = stats(window('burn_start', 'burn_end', 5.0)); idle2 = stats(window('burn_end', 'idle_end'))
    bj = json.load(open('results/W/burn.burn.json'))
    out['W'] = dict(idle=idle, burn=burn, idle_after=idle2, burn_run=bj, clock_offset_s=off.total_seconds())
    md.append('\n## Power (GPU1 by UUID, 200 ms sampling)\n\n| Window | Samples | util % | board W (power.draw) | board W (instant) | module W (avg) | module W (instant) | SM MHz |\n|---|---:|---:|---:|---:|---:|---:|---:|')
    for k, s in [('idle before', idle), ('burn (75 s back-to-back batches)', burn), ('idle after', idle2)]:
        md.append(f"| {k} | {s['samples']} | {fmt(s['util'], 1)} | {fmt(s['power_draw_w'], 1)} | {fmt(s['power_draw_instant_w'], 1)} | {fmt(s['module_power_avg_w'], 1)} | {fmt(s['module_power_instant_w'], 1)} | {fmt(s['sm_clock_mhz'], 0)} |")
    md.append(f"\nBurn run: {bj['batches']} batches of B = {bj['chains']} in {bj['burn_seconds_wall']:.1f} s, {bj['device_ms_per_batch']:.4f} ms per batch.")
json.dump(out, open('results/tables.json', 'w'), indent=1, default=float)
open('results/tables.md', 'w').write('\n'.join(md) + '\n')
print('\n'.join(md))
