"""Analysis of PROTOCOL.md results -> results_summary.json, selected_configs.json, analysis.txt (tables used in REPORT.md).
Reads results/budget/G*/<method>_S<S>.json and results/conf/G*/<rule>_<E1|E12>.json. Usage: python3 analyze_gset.py"""
import json
import math
import sys
from pathlib import Path

import numpy as np

sys.dont_write_bytecode = True
import engine as E  # noqa: E402
import run_gset as R  # noqa: E402

ROOT = E.ROOT
RES = ROOT / 'results'
MAN = json.loads((ROOT / 'manifest.json').read_text())['instances']
BKV = R.BKV
CLASS = {g: MAN[f'G{g}']['cls'] for g in R.ALL}
CLASS_ORDER = ['R800+', 'R800+-', 'T800+-', 'P800+', 'P800+-', 'R2000+', 'R2000+-', 'T2000+-', 'P2000+', 'P2000+-', 'R1000+', 'P1000+']
RULES = R.ENGINE_RULES
SHORT = {'SCA': 'plain', 'TEC': 'TEC', 'Onsager-kT': 'kT', 'Onsager-online': 'online'}
NBOOT = 2000


VARIANT = 'orig'      # 'orig' = frozen PROTOCOL.md grids; 'ext' = PROTOCOL_AMENDMENT1.md union grids (engine rules only)


def ext_budget_path(g, m, S):
    return RES / 'ext_budget' / f'G{g}' / f"{m.replace(' ', '_')}_S{S}.json"


def ext_conf_path(g, m, o):
    return RES / 'ext_conf' / f'G{g}' / f"{m.replace(' ', '_')}_{o}.json"


def budget_file(g, m, S):
    return ext_budget_path(g, m, S) if (VARIANT == 'ext' and m in RULES) else R.budget_path(g, m, S)


def conf_file(g, m, o):
    return ext_conf_path(g, m, o) if VARIANT == 'ext' else R.conf_path(g, m, o)


def load_budget(g, m, S):
    p = budget_file(g, m, S)
    if not p.exists():
        return None
    d = json.loads(p.read_text())
    if 'union_selected_cfg' in d:
        d['selected_cfg'] = d['union_selected_cfg']; d['selected'] = d['union_selected']
    return d


def load_conf(g, m, o):
    p = conf_file(g, m, o)
    return json.loads(p.read_text()) if p.exists() else None


def tts_from_runs(cuts, flips, S, target, obj):
    p = float((cuts >= target).mean())
    if obj == 'E1':
        return E.t_trial_ms(float(np.mean(flips)), S) * E.r99(p)
    tr = E.expected_max(E.t_trial_ms(flips, S), E.E_ENGINES); P = 1 - (1 - p) ** E.E_ENGINES
    return tr if P >= 0.995 else tr * E.r99(P)


def conf_tts(g, m, o):
    d = load_conf(g, m, o)
    if d is None:
        return None, None
    if d['final'] is None:
        return math.inf, d
    return d['final']['tts1_ms' if o == 'E1' else 'tts12_ms'], d


def boot_ratio(dref, dx, o, g, rng):
    """Paired bootstrap (same run indices; confirmation seeds are shared across rules) of TTS_ref / TTS_x."""
    tgt = BKV[f'G{g}']['target']
    if dref is None or dx is None or dref['final'] is None or dx['final'] is None:
        return None
    cr, fr = np.array(dref['cuts']), np.array(dref['flips']); cx, fx = np.array(dx['cuts']), np.array(dx['flips'])
    n = len(cr); out = np.empty(NBOOT)
    for b in range(NBOOT):
        idx = rng.integers(0, n, n)
        a = tts_from_runs(cr[idx], fr[idx], dref['S'], tgt, o); x = tts_from_runs(cx[idx], fx[idx], dx['S'], tgt, o)
        if math.isfinite(x):
            out[b] = a / x                      # a may be inf -> inf (only the rule solves in this resample)
        elif math.isfinite(a):
            out[b] = 0.0                        # only the reference solves
        else:
            out[b] = np.nan                     # neither solves: undefined
    with np.errstate(invalid='ignore'):
        return [float(np.nanpercentile(out, 2.5)), float(np.nanpercentile(out, 97.5))]


def ratio(a, x):
    """speedup of x over reference a (TTS_a / TTS_x); inf if only x solves, 0 if only a solves, None if neither."""
    if a is None or x is None:
        return None
    if math.isinf(a) and math.isinf(x):
        return None
    if math.isinf(x):
        return 0.0
    if math.isinf(a):
        return math.inf
    return a / x


def gstats(vals):
    v = [x for x in vals if x is not None and math.isfinite(x) and x > 0]
    if not v:
        return None
    lv = np.log(v)
    return dict(n=len(v), geomean=float(np.exp(lv.mean())), min=float(min(v)), median=float(np.median(v)), max=float(max(v)),
                frac_gt1=float(np.mean(np.array(v) > 1)))


def fmt(x, d=3):
    if x is None:
        return '–'
    if isinstance(x, float) and math.isinf(x):
        return '∞'
    return f'{x:.{d}g}'


def main():
    rng = np.random.default_rng(777)
    inst = {}
    have = [g for g in R.ALL if all(budget_file(g, m, S).exists() for m in RULES for S in R.S_LIST)]
    conf_done = [g for g in have if all(conf_file(g, m, o).exists() for m in RULES for o in R.OBJECTIVES)]
    for g in have:
        b = BKV[f'G{g}']
        rec = dict(cls=CLASS[g], N=MAN[f'G{g}']['N'], m=MAN[f'G{g}']['m'], sigma=MAN[f'G{g}']['sigma'], BKV=b['BKV'], target=b['target'],
                   rules={})
        for m in R.METHODS:
            rows = {S: load_budget(g, m, S) for S in R.S_LIST}
            if not all(rows.values()):
                continue
            per_S = {}
            for S, d in rows.items():
                f = d['final']
                per_S[S] = dict(selected=d['selected'], cfg=d['selected_cfg'], **({} if f is None else
                                {k: f[k] for k in f if k in ('mean_cut', 'sd_cut', 'max_cut', 'p_target', 'p_target_wilson95', 'k_target',
                                                             'k_bkv', 'runs', 'mcs99', 'tts1_ms', 'tts12_ms', 'mean_flips',
                                                             't_trial_ms', 't_round12_ms', 'n_invalid', 'k_target_13221')}))
            mcs = [(per_S[S].get('mcs99', math.inf), S) for S in R.S_LIST]
            best = min(mcs)
            rr = dict(per_S=per_S, mcs99_min=best[0], mcs99_S=best[1] if math.isfinite(best[0]) else None,
                      k_bkv_finals=sum(per_S[S].get('k_bkv', 0) for S in R.S_LIST),
                      runs_finals=sum(per_S[S].get('runs', 0) for S in R.S_LIST),
                      best_cut_finals=max((per_S[S].get('max_cut') or -1) for S in R.S_LIST))
            if m in RULES and g in conf_done:
                for o in R.OBJECTIVES:
                    v, d = conf_tts(g, m, o)
                    rr[f'conf_{o}'] = None if d is None else dict(S=d['S'], cfg=d['cfg'], tts_ms=v,
                                                                  **({} if d['final'] is None else {k: d['final'][k] for k in (
                                                                      'p_target', 'p_target_wilson95', 'k_target', 'k_bkv', 'runs',
                                                                      'mean_cut', 'max_cut', 'mean_flips', 't_trial_ms',
                                                                      't_round12_ms', 'P_round12') + (('k_target_13221',) if g == 23 else ())}))
            rec['rules'][m] = rr
        if g in conf_done:
            rec['ratios'] = {}
            for o in R.OBJECTIVES:
                t = {m: conf_tts(g, m, o) for m in RULES}
                for ref in ('SCA', 'TEC'):
                    for x in ('Onsager-kT', 'Onsager-online') + (('TEC',) if ref == 'SCA' else ()):
                        rec['ratios'][f'{o}:{SHORT[ref]}/{SHORT[x]}'] = dict(
                            ratio=ratio(t[ref][0], t[x][0]), ci95=boot_ratio(t[ref][1], t[x][1], o, g, rng))
                bo = min(t['Onsager-kT'][0], t['Onsager-online'][0])
                rec['ratios'][f'{o}:plain/bestOnsager'] = dict(ratio=ratio(t['SCA'][0], bo))
                rec['ratios'][f'{o}:TEC/bestOnsager'] = dict(ratio=ratio(t['TEC'][0], bo))
        inst[f'G{g}'] = rec
    # hypotheses (only meaningful when all 51 instances are confirmed)
    hyp = {}
    for o in R.OBJECTIVES:
        for hid, ref, x in (('G1', 'SCA', 'Onsager-online'), ('G2', 'SCA', 'Onsager-kT'), ('G3', 'TEC', 'Onsager-online'),
                            ('G4', 'TEC', 'Onsager-kT')):
            wins = []; ties = []; losses = []
            for gname, rec in inst.items():
                if 'ratios' not in rec:
                    continue
                a = rec['rules'][ref][f'conf_{o}']['tts_ms']; b = rec['rules'][x][f'conf_{o}']['tts_ms']
                (wins if b < a else ties if b == a else losses).append(gname)
            n = len(wins) + len(ties) + len(losses)
            hyp[f'{hid}{"" if o == "E1" else chr(39)}'] = dict(objective=o, reference=ref, rule=x, wins=len(wins), ties=len(ties),
                                                               losses=len(losses), instances=n,
                                                               pass_=bool(n == 51 and len(wins) > n / 2), win_list=wins)
        gm = gstats([rec['ratios'][f'{o}:plain/bestOnsager']['ratio'] for rec in inst.values() if 'ratios' in rec])
        hyp[f'G5_{o}'] = dict(objective=o, geomean_plain_over_bestOnsager=gm, pass_=bool(gm is not None and gm['geomean'] >= 2
                                                                                          and len(conf_done) == 51))
    summ = {}
    for o in R.OBJECTIVES:
        for key in (f'{o}:plain/online', f'{o}:plain/kT', f'{o}:plain/TEC', f'{o}:TEC/online', f'{o}:TEC/kT',
                    f'{o}:plain/bestOnsager', f'{o}:TEC/bestOnsager'):
            vals = [rec['ratios'][key]['ratio'] for rec in inst.values() if 'ratios' in rec]
            summ[key] = dict(stats=gstats(vals), n_x_only=sum(1 for v in vals if v == math.inf),
                             n_ref_only=sum(1 for v in vals if v == 0.0), n_neither=sum(1 for v in vals if v is None))
            percls = {}
            for c in CLASS_ORDER:
                vv = [rec['ratios'][key]['ratio'] for rec in inst.values() if 'ratios' in rec and rec['cls'] == c]
                if vv:
                    percls[c] = dict(n=len(vv), x_wins=sum(1 for v in vv if v is not None and v > 1), stats=gstats(vv))
            summ[key]['per_class'] = percls
    # steps to solution (hardware-independent): MCS99 min over S on the 256-run finals
    mcs = {}
    for ref, x in (('SCA', 'Onsager-online'), ('SCA', 'Onsager-kT'), ('SCA', 'TEC'), ('TEC', 'Onsager-online'), ('TEC', 'Onsager-kT')):
        vals = [ratio(rec['rules'][ref]['mcs99_min'], rec['rules'][x]['mcs99_min']) for rec in inst.values()
                if ref in rec['rules'] and x in rec['rules']]
        percls = {}
        for c in CLASS_ORDER:
            vv = [ratio(rec['rules'][ref]['mcs99_min'], rec['rules'][x]['mcs99_min']) for rec in inst.values()
                  if rec['cls'] == c and ref in rec['rules'] and x in rec['rules']]
            if vv:
                percls[c] = dict(n=len(vv), x_lower=sum(1 for v in vv if v is not None and v > 1), stats=gstats(vv))
        mcs[f'{SHORT[ref]}/{SHORT[x]}'] = dict(stats=gstats(vals), n_x_only=sum(1 for v in vals if v == math.inf),
                                               n_ref_only=sum(1 for v in vals if v == 0.0), n_neither=sum(1 for v in vals if v is None),
                                               x_lower=sum(1 for v in vals if v is not None and v > 1),
                                               x_higher=sum(1 for v in vals if v is not None and v < 1), n=len(vals), per_class=percls)
    sub = [f'G{g}' for g in R.SUBSET if f'G{g}' in inst and all(m in inst[f'G{g}']['rules'] for m in R.METHODS)]
    allm = {}
    for m in R.METHODS:
        v = [inst[g]['rules'][m]['mcs99_min'] for g in sub]
        rel = [ratio(inst[g]['rules'][m]['mcs99_min'], inst[g]['rules']['Onsager-online']['mcs99_min']) for g in sub]
        allm[m] = dict(solved=sum(1 for x in v if math.isfinite(x)), of=len(sub),
                       geomean_mcs_over_online=gstats(rel),
                       bkv_hits=sum(inst[g]['rules'][m]['k_bkv_finals'] for g in sub),
                       bkv_instances=sum(1 for g in sub if inst[g]['rules'][m]['k_bkv_finals'] > 0),
                       runs=sum(inst[g]['rules'][m]['runs_finals'] for g in sub))
    out = dict(variant=VARIANT, instances_with_partA=len(have), instances_confirmed=len(conf_done), hypotheses=hyp, summary=summ,
               mcs99_summary_engine_rules=mcs, subset_all_methods=dict(instances=sub, per_method=allm), instances=inst)
    suf = '' if VARIANT == 'orig' else '_ext'
    (ROOT / f'results_summary{suf}.json').write_text(json.dumps(out, indent=1, default=lambda x: None) + '\n')
    write_text(out, suf)
    return inst


def write_selected(inst_orig, inst_ext):
    sel = dict(note=('Per instance, the selected configuration of each SCA-family rule at every budget S (pilot mean-cut rule), and the '
                     'configurations confirmed for single-engine (E1) and 12-engine (E12) TTS. Absolute values are in the units of the '
                     'G-set couplings J = -w (weights +-1). Rule: z_i = s_i*field_i + q; stay probability clip(z/(4T) + 1/2, 0, 1); '
                     'T(t) = T0*(tfin/T0)^(t/(S-1)); plain: field = h; TEC: field = h + jv*s(t-1); Onsager-kT (tecT): field = h - '
                     'kappa*T(t-1)*ramp(t)*s(t-1); Onsager-online (onsager): field = h - lam*ramp(t)*c(t-1)*s(t-1), c = n_lin/(2T), '
                     'n_lin = #{i: -2T < z_i < 2T} over the N active spins (exclude padding slots); lam already includes the factor '
                     '(dbar/N)/(1999/2000) (rel.lam_k2000 = lambda in K2000 units); APC: per-spin q_i reset to q_reset after a flip, '
                     'else q_i <- max(r_q*q_i, q_lim). ramp(t) = 1 for t < 0.7S, then (S-1-t)/max(1, S-1-0.7S). Model: '
                     'research/ablation_20261005/abl.py::run (engine.py is bit-identical to it).'), instances={})
    for gname, rec in inst_orig.items():
        e = {}
        for variant, inst in (('original', inst_orig), ('extended_amendment1', inst_ext)):
            if inst is None or gname not in inst:
                continue
            for m in R.ENGINE_RULES + ('APC-SCA',):
                if m not in inst[gname]['rules'] or (variant != 'original' and m == 'APC-SCA'):
                    continue
                rr = inst[gname]['rules'][m]
                v = e.setdefault(m, {}).setdefault(variant, dict(per_budget={str(S): rr['per_S'][S]['cfg'] for S in R.S_LIST}))
                for o in R.OBJECTIVES:
                    if f'conf_{o}' in rr and rr[f'conf_{o}']:
                        v[f'tts_{o}'] = dict(S=rr[f'conf_{o}']['S'], cfg=rr[f'conf_{o}']['cfg'], tts_ms=rr[f'conf_{o}']['tts_ms'],
                                             p=rr[f'conf_{o}'].get('p_target'))
        sel['instances'][gname] = dict(N=rec['N'], target=rec['target'], BKV=rec['BKV'], rules=e)
    (ROOT / 'selected_configs.json').write_text(json.dumps(sel, indent=1, default=lambda x: None) + '\n')


def write_text(out, suf=''):
    L = []
    inst = out['instances']
    L.append(f"variant {out['variant']}: instances with Part A complete: {out['instances_with_partA']}; confirmed (Part C): {out['instances_confirmed']}")
    L.append('\nPer-instance engine TTS (confirmation runs, 1024 each). p = P(cut >= target); TTS in ms; ratio = TTS_plain/TTS_rule')
    hdr = (f"{'inst':5s} {'class':8s} {'BKV':>6s} {'tgt':>6s} | " + ' | '.join(f"{SHORT[m]:>6s} S  p     TTS1" for m in RULES)
           + ' | plain/kT plain/onl TEC/kT TEC/onl (E1) | same, E12')
    L.append(hdr)
    for gname, rec in inst.items():
        if 'ratios' not in rec:
            continue
        cells = []
        for m in RULES:
            c = rec['rules'][m]['conf_E1']
            cells.append(f"{fmt(c['S'], 4):>6s} {fmt(c.get('p_target'), 2):>4s} {fmt(c['tts_ms'], 3):>7s}")
        r = rec['ratios']
        rs = ' '.join(f"{fmt(r[k]['ratio'], 3):>7s}" for k in ('E1:plain/kT', 'E1:plain/online', 'E1:TEC/kT', 'E1:TEC/online'))
        r12 = ' '.join(f"{fmt(r[k]['ratio'], 3):>7s}" for k in ('E12:plain/kT', 'E12:plain/online', 'E12:TEC/kT', 'E12:TEC/online'))
        L.append(f"{gname:5s} {rec['cls']:8s} {rec['BKV']:6d} {rec['target']:6d} | " + ' | '.join(cells) + f" | {rs} | {r12}")
    L.append('\nHypotheses')
    for k, v in out['hypotheses'].items():
        if k.startswith('G5'):
            L.append(f"  {k}: geomean plain/best-Onsager = {v['geomean_plain_over_bestOnsager']} pass={v['pass_']}")
        else:
            L.append(f"  {k} ({v['objective']}): {v['rule']} lower TTS than {v['reference']} on {v['wins']}/{v['instances']} "
                     f"(ties {v['ties']}, losses {v['losses']}) pass={v['pass_']}")
    L.append('\nSummary of TTS ratios (geometric mean over instances where both finite; n_x_only = only the rule solves)')
    for k, v in out['summary'].items():
        s = v['stats']
        L.append(f"  {k:22s} " + (f"geomean {s['geomean']:.3f} [min {s['min']:.3f}, median {s['median']:.3f}, max {s['max']:.3f}] "
                                  f"n={s['n']} frac>1 {s['frac_gt1']:.2f}" if s else 'n/a') +
                 f" | x-only {v['n_x_only']} ref-only {v['n_ref_only']} neither {v['n_neither']}")
    L.append('\nMCS99 ratios (ref/rule, min over S) for the engine rules on all Part-A instances')
    for k, v in out['mcs99_summary_engine_rules'].items():
        s = v['stats']
        L.append(f"  {k:14s} " + (f"geomean {s['geomean']:.3f} [min {s['min']:.3f}, median {s['median']:.3f}, max {s['max']:.3f}] n={s['n']}"
                                  if s else 'n/a') + f" | rule lower {v['x_lower']} higher {v['x_higher']} of {v['n']}; "
                                                     f"rule-only {v['n_x_only']} ref-only {v['n_ref_only']} neither {v['n_neither']}")
    sa = out['subset_all_methods']
    L.append(f"\nAll ten methods on the subset ({len(sa['instances'])} instances complete): solved, geomean MCS99 / Onsager-online, BKV hits")
    for m, v in sa['per_method'].items():
        s = v['geomean_mcs_over_online']
        L.append(f"  {m:15s} solved {v['solved']}/{v['of']} | MCS99 ratio to online: " +
                 (f"geomean {s['geomean']:.3f} [min {s['min']:.3g}, max {s['max']:.3g}] n={s['n']}" if s else 'n/a') +
                 f" | BKV hits {v['bkv_hits']}/{v['runs']} runs on {v['bkv_instances']} instances")
    L.append('\nSteps to solution MCS99 (min over S, 256-run finals) and BKV hits (all finals)')
    meths = [m for m in R.METHODS]
    L.append(f"{'inst':5s} " + ' '.join(f"{m[:10]:>11s}" for m in meths))
    for gname, rec in inst.items():
        L.append(f"{gname:5s} " + ' '.join(f"{fmt(rec['rules'][m]['mcs99_min'], 4) if m in rec['rules'] else '':>11s}" for m in meths))
    L.append('\nBKV hits / runs in the 5 finals')
    for gname, rec in inst.items():
        L.append(f"{gname:5s} " + ' '.join(f"{(str(rec['rules'][m]['k_bkv_finals']) + '/' + str(rec['rules'][m]['runs_finals'])) if m in rec['rules'] else '':>11s}"
                                         for m in meths))
    (ROOT / f'analysis{suf}.txt').write_text('\n'.join(L) + '\n')
    print('\n'.join(L))


if __name__ == '__main__':
    inst_o = main()
    inst_e = None
    if (RES / 'ext_budget').exists():
        VARIANT = 'ext'
        inst_e = main()
    write_selected(inst_o, inst_e)
