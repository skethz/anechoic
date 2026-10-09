"""Print Markdown tables built from results/*.json (stdout; saved as logs/summarize.log). Usage: python3 summarize.py"""
import json
import math
from pathlib import Path

ROOT = Path(__file__).resolve().parent
RES = ROOT / 'results'
GBSB = 9.61
FR = ('V80_250MHz', 'ASIC_1.0GHz', 'ASIC_1.25GHz')


def ld(n):
    p = RES / f'{n}.json'
    return json.loads(p.read_text()) if p.exists() else None


def lab(c):
    f = {'onsager': 'Onsager-online', 'tecT': 'Onsager-kT', 'plain': 'Plain'}[c['family']]
    par = f"q{c['q']:g}"
    if 'lam' in c: par += f" λ{c['lam']:g}"
    if 'kappa' in c: par += f" κ{c['kappa']:g}"
    s = f"{f} S{c['S']} T0={c['T0']:g} T1={c.get('T1', 5):g} {par}"
    if c.get('fin_L'):
        s += f" + fin L{c['fin_L']} (q{c['fin_q']:g}, T{c['fin_T']:g})"
    return s


def ms(x):
    return '∞' if x is None or not math.isfinite(x) else (f'{x:.1f}' if x >= 10 else f'{x:.2f}')


def pci(k, n, lo, hi):
    return f"{k}/{n} = {k / n:.4f} [{lo:.4f}, {hi:.4f}]"


def wilson(k, n, z=1.959963984540054):
    p = k / n; d = 1 + z * z / n
    c = (p + z * z / (2 * n)) / d; w = z * math.sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / d
    return max(0.0, c - w), min(1.0, c + w)


def tts_cell(t, fr):
    x = t[fr]
    return f"{ms(x['tts12_ms'])} [{ms(x['tts12_ms_lo'])}, {ms(x['tts12_ms_hi'])}]"


def main():
    out = []
    a1 = ld('a1_screen')
    if a1:
        rows = a1['rows']
        out.append('## (a) Stage 1 screen (256 runs per configuration)\n')
        out.append('| Family | S | configs | best mean final cut (config) | max k(≥33,300) | Σ k(33,337) | Σ k_greedy(33,337) |')
        out.append('|---|---:|---:|---|---:|---:|---:|')
        groups = {}
        for r in rows:
            groups.setdefault((r['cfg']['family'], r['cfg']['S']), []).append(r)
        for (fam, S), rs in sorted(groups.items(), key=lambda x: (['onsager', 'tecT', 'plain'].index(x[0][0]), x[0][1])):
            b = max(rs, key=lambda r: r['mean_cut'])
            out.append(f"| {fam} | {S} | {len(rs)} | {b['mean_cut']:.1f} ({lab(b['cfg'])}) | {max(r['k33300'] for r in rs)} | "
                       f"{sum(r['k33337'] for r in rs)} | {sum(r['kgreedy33337'] for r in rs)} |")
        # T1 effect: paired by all other parameters
        out.append('\n**Effect of the final temperature T1 inside the geometric schedule** (mean final cut, averaged over all '
                   'other grid parameters; corrected families exclude runaway configurations with mean cut < 33,050):\n')
        out.append('| Family | S | T1 = 5 | T1 = (q+1)/2 | T1 = q/2 |')
        out.append('|---|---:|---:|---:|---:|')
        for (fam, S), rs in sorted(groups.items(), key=lambda x: (['onsager', 'tecT', 'plain'].index(x[0][0]), x[0][1])):
            base = {}
            for r in rs:
                c = dict(r['cfg']); t1 = c.pop('T1'); q = c['q']
                lvl = 0 if t1 == 5.0 else (1 if t1 == (q + 1) / 2 else 2)
                base.setdefault(json.dumps(c, sort_keys=True), {})[lvl] = r['mean_cut']
            good = [v for v in base.values() if len(v) == 3 and (fam == 'plain' or v[0] >= 33050)]
            if good:
                m = [sum(v[l] for v in good) / len(good) for l in range(3)]
                out.append(f"| {fam} | {S} | {m[0]:.1f} | {m[1]:.1f} | {m[2]:.1f} |")
        out.append('\n**Final-state vs host-greedy vs best-visited hits at 33,337, summed over the screen** (256 runs per '
                   'configuration; descriptive):\n')
        out.append('| Family | end point | configs | runs | final k | greedy k | best-visited k | final/best-visited |')
        out.append('|---|---|---:|---:|---:|---:|---:|---:|')
        for fam in ('onsager', 'tecT', 'plain'):
            for lvl, name in ((0, 'T1 = 5'), (1, 'T1 = (q+1)/2'), (2, 'T1 = q/2')):
                rs = [r for r in rows if r['cfg']['family'] == fam and
                      (0 if r['cfg']['T1'] == 5.0 else (1 if r['cfg']['T1'] == (r['cfg']['q'] + 1) / 2 else 2)) == lvl]
                k = sum(r['k33337'] for r in rs); kg = sum(r['kgreedy33337'] for r in rs); kb = sum(r['kbest33337'] for r in rs)
                out.append(f"| {fam} | {name} | {len(rs)} | {sum(r['n'] for r in rs)} | {k} | {kg} | {kb} | "
                           f"{(k / kb) if kb else float('nan'):.2f} |")
    a2 = ld('a2_confirm')
    if a2:
        out.append('\n## (a) Stage 2 confirm-pilot (2,048 runs per candidate)\n')
        out.append('| Candidate | mean cut | k≥33,300 | k≥33,320 | k=33,337 | k_greedy 33,337 | k_best-visited 33,337 | t_trial V80 [ms] | score [ms] |')
        out.append('|---|---:|---:|---:|---:|---:|---:|---:|---:|')
        for r in a2['rows']:
            lo, _ = wilson(r['k33337'], r['n'])
            tR = r['round_cycles_E12'] / 250e6
            sc = 1e3 * tR * math.log(0.01) / (12 * math.log1p(-lo)) if lo > 0 else math.inf
            out.append(f"| {lab(r['cfg'])} | {r['mean_cut']:.1f} | {r['k33300']} | {r['k33320']} | {r['k33337']} | {r['kgreedy33337']} | "
                       f"{r['kbest33337']} | {r['tts']['V80_250MHz']['t_trial_ms']:.3f} | {ms(sc)} |")
        if 'selection' in a2:
            out.append('\nSelection: ' + '; '.join(f"**{k}**: {lab(v['cfg'])} (ranked at {v['ranked_at']}, score {ms(v['score_ms'])} ms)"
                                                  for k, v in a2['selection'].items()))
    cp = ld('c_pilot')
    if cp:
        out.append('\n## (c) Finishing pilot (4,096 anneals per configuration; same anneals for every option)\n')
        for kk, v in cp['configs'].items():
            rows = v['rows']; none = rows[0]
            out.append(f"\n**{lab(v['cfg'])}** — none: k337 {none['k33337']}, greedy {none['kgreedy33337']}, best-visited "
                       f"{none['kbest33337']}; selected {v['selected']} (score {ms(v['selected_score_ms'])} vs none {ms(v['none_score_ms'])} ms)\n")
            out.append('| L | (q_f, T_f) | k≥33,300 | k≥33,320 | k=33,337 | k_greedy 33,337 | Δ cycles/trial |')
            out.append('|---:|---|---:|---:|---:|---:|---:|')
            for r in rows:
                o = r['opt']
                d = r['trial_cycles'] - none['trial_cycles']
                out.append(f"| {o['fin_L'] if o else 0} | {('(%g, %g)' % (o['fin_q'], o['fin_T'])) if o else '—'} | {r['k33300']} | "
                           f"{r['k33320']} | {r['k33337']} | {r['kgreedy33337']} | {d:.0f} |")
    h = ld('holdout_ac')
    if h:
        out.append('\n## Held-out (8,192 runs per configuration, paired starts): final-state success at 33,337\n')
        out.append('| Configuration | S_total | p(33,337) [95% CI] | best-visited k | greedy k | t_R (12 eng., V80) [ms] | '
                   'TTS99 V80 12 eng. [ms] | ASIC 1.0 GHz | ASIC 1.25 GHz | V80 / GbSB |')
        out.append('|---|---:|---|---:|---:|---:|---|---|---|---:|')
        for kk, v in h['configs'].items():
            for r in [v['no_finish']] + ([v['selected']] if 'selected' in v else []):
                t = r['tts']
                out.append(f"| {lab(r['cfg'])} | {r['S_total']} | {pci(r['k33337'], r['n'], t['p_lo'], t['p_hi'])} | {r['kbest33337']} | "
                           f"{r['kgreedy33337']} | {t['V80_250MHz']['t_round_ms']:.3f} | {tts_cell(t, 'V80_250MHz')} | "
                           f"{tts_cell(t, 'ASIC_1.0GHz')} | {tts_cell(t, 'ASIC_1.25GHz')} | {t['V80_250MHz']['tts12_ms'] / GBSB:.1f}× |")
        out.append('\n**Paired finishing comparison (selected option, same 8,192 anneals):**\n')
        out.append('| Configuration | none k | finish k | both | finish only | anneal only | added cycles/trial | Δ TTS V80 |')
        out.append('|---|---:|---:|---:|---:|---:|---:|---|')
        for kk, v in h['configs'].items():
            if 'selected' not in v:
                continue
            a, b = v['no_finish'], v['selected']
            pr = v['paired']
            out.append(f"| {lab(a['cfg'])} + {v['selected_option']} | {a['k33337']} | {b['k33337']} | {pr['both']} | {pr['finish_only']} | "
                       f"{pr['anneal_only']} | {b['trial_cycles'] - a['trial_cycles']:.0f} ({100 * (b['trial_cycles'] / a['trial_cycles'] - 1):.1f}%) | "
                       f"{ms(a['tts']['V80_250MHz']['tts12_ms'])} → {ms(b['tts']['V80_250MHz']['tts12_ms'])} ms |")
        out.append('\n**All finishing options on the held-out anneals (exploratory except the selected one):**\n')
        for kk, v in h['configs'].items():
            a = v['no_finish']
            best = sorted(v['options'], key=lambda r: r['tts']['V80_250MHz']['tts12_ms'])[:5]
            out.append(f"- {lab(a['cfg'])}: none k={a['k33337']} (greedy {a['kgreedy33337']}); top-5 options by TTS: " +
                       '; '.join(f"L{r['opt']['fin_L']} (q{r['opt']['fin_q']:g},T{r['opt']['fin_T']:g}) k={r['k33337']} "
                                 f"TTS {ms(r['tts']['V80_250MHz']['tts12_ms'])} ms" for r in best))
    bp = ld('b_pilot')
    if bp:
        out.append('\n## (b) Population annealing pilot (256 populations of R = 12)\n')
        out.append(f"Base: {json.dumps(bp['base'])}\n")
        out.append('| S | variant | K | α | P_pop(33,337) | P_pop(33,320) | P_pop(33,300) | single-replica k337 | copied/event | t_pop V80 [ms] | TTS V80 [ms] |')
        out.append('|---:|---|---:|---:|---:|---:|---:|---:|---:|---:|---|')
        for r in bp['rows']:
            c = r['cfg']; t = r['tts']['main|V80_250MHz']
            out.append(f"| {c['S']} | {c['pa']} | {c['K']} | {c['alpha']:g} | {r['k33337']}/{r['populations']} | {r['k33320']} | {r['k33300']} | "
                       f"{r['single_replica_k33337']}/{r['single_replica_n']} | {r['mean_copied_per_event']:.2f} | {t['t_pop_ms']:.3f} | {ms(t['tts_ms'])} |")
        if 'selection' in bp:
            out.append(f"\nSelected: {json.dumps(bp['selection'])}")
    for name in ('holdout_b', 'holdout_b48'):
        hb = ld(name)
        if not hb:
            continue
        out.append(f"\n## (b) Held-out {name}\n")
        out.append('| Variant | R | populations | P_pop(33,337) [95% CI] | single-replica p | P_pop greedy | P_pop best-visited | cost variant | t_pop [ms] | TTS V80 | ASIC 1.0 GHz | ASIC 1.25 GHz |')
        out.append('|---|---:|---:|---|---:|---:|---:|---|---:|---|---|---|')
        for v, r in hb['results'].items():
            for cv in r['pop_cycles']:
                for tag in ('', '|R_engines'):
                    keyv = f'{cv}|V80_250MHz{tag}'
                    if keyv not in r['tts']:
                        continue
                    t = {fr: r['tts'][f'{cv}|{fr}{tag}'] for fr in FR}
                    cell = lambda fr: f"{ms(t[fr]['tts_ms'])} [{ms(t[fr]['tts_ms_lo'])}, {ms(t[fr]['tts_ms_hi'])}]"
                    out.append(f"| {v} | {r['R']} | {r['populations']} | {pci(r['k33337'], r['populations'], r['P_pop_lo'], r['P_pop_hi'])} | "
                               f"{r['single_replica_k33337'] / r['single_replica_n']:.4f} | {r['kgreedy33337']} | {r['kbest33337']} | "
                               f"{cv}{' (R engines)' if tag else (' (12 engines)' if r['R'] > 12 else '')} | {t['V80_250MHz']['t_pop_ms']:.3f} | "
                               f"{cell('V80_250MHz')} | {cell('ASIC_1.0GHz')} | {cell('ASIC_1.25GHz')} |")
    for name in ('explore_e1', 'explore_e2'):
        e = ld(name)
        if e:
            out.append(f"\n## Exploratory {name}\n")
            out.append('```\n' + json.dumps({k: {kk: {x: vv[x] for x in ('n', 'mean_cut', 'k33300', 'k33320', 'k33337', 'kgreedy33337', 'S_total')
                                                       if x in vv} | {'tts12_V80_ms': vv['tts']['V80_250MHz']['tts12_ms']}
                                                  for kk, vv in v.items()} for k, v in e['results'].items()}, indent=1) + '\n```')
    print('\n'.join(out))


if __name__ == '__main__':
    main()
