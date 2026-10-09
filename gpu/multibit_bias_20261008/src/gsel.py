#!/usr/bin/env python3
"""PROTOCOL.md Amendment 5: GPU-selected G-set schedules, the GPU counterpart of the study's V80 selection.

The study chose the V80's E12 configuration as argmin over S in {250, 500, 1000, 2000, 4000} of its V80 12-engine round
estimator evaluated on the 256-run finals of the per-budget (mean-cut) selected configuration (verified for all 408
original and extended E12/E1 choices). This script applies the same procedure with a GPU cost model:

  candidates   (S, B), S in the five budgets, B in {60, 120, 240}; the configuration at S is the study's per-budget
               selection (original grid: budget/<g>/<rule>_S<S>.json selected_cfg; extended grid: ext_budget/... union
               selection), with its 256-run final (k successes of n = 256)
  cost model   t(S, B, class) = median measured device time per batch of the phase-G G-set runs with that (B, class, S)
               (class = const for SCA, TEC, Onsager-kT; ons for Onsager-online), i.e. the measured batch-time floor of
               the variant the driver uses for that (B, class, S)
  primary      P_B = 1 - (1 - p)^B; TTS = t if P_B >= 0.995 else t ln(0.01) / ln(1 - P_B)  (the study's round convention)
  secondary    TTS = t ln(0.01) / (B ln(1 - p))
               p = k / n, with k = n replaced by n - 0.5 (continuity); k = 0 gives TTS = inf
  selection    argmin of the primary prediction, and separately argmin of the secondary prediction, ties to smaller S
               then smaller B; per instance, rule (SCA, TEC, Onsager-kT, Onsager-online) and grid (original, extended)
Also computes the V80's own secondary-optimal selection with the study's V80 model (argmin over S of
t_round12 ln(0.01)/(12 ln(1-p)) on the same finals), used only for a labelled modelled comparison.
Writes results/GS/plan_labels.json, results/GS/plan_cohorts.tsv and results/GS/cost_model.json. Deterministic."""
import glob, json, math, os, statistics, sys

RULES = ['SCA', 'TEC', 'Onsager-kT', 'Onsager-online']
CLS = {'SCA': 'const', 'TEC': 'const', 'Onsager-kT': 'const', 'Onsager-online': 'ons'}
SS = [250, 500, 1000, 2000, 4000]
BS = [60, 120, 240]
LN01 = math.log(0.01)
STUDY = sys.argv[1] if len(sys.argv) > 1 else 'data/gset_study_results'


def flags(cfg):   # identical mapping to gset_plan.py
    f = cfg['family']
    common = f"--q {cfg['q']!r} --t0 {cfg['T0']!r} --t1 {cfg['tfin']!r} --steps {int(cfg['S'])}"
    if f == 'plain': return '--lambda 0 ' + common
    if f == 'tec': return f"--tec-jv {cfg['jv']!r} " + common
    if f == 'tecT': assert cfg.get('ramp') is True; return f"--tecT {cfg['kappa']!r} --ramp " + common
    if f == 'onsager': assert cfg.get('ramp') is True; return f"--lambda {cfg['lam']!r} --ramp " + common
    raise ValueError(f)


def cost_model():
    plan = {c['tag']: c for c in json.load(open('results/G/plan.json'))}
    grp = {}
    for f in sorted(glob.glob('results/G/*_b*.summary.json')):
        s = json.load(open(f)); B = s['chains']
        if B == 1: continue
        c = plan[os.path.basename(f)[:-len('.summary.json')].rsplit('_b', 1)[0]]
        grp.setdefault(f"{B}|{c['cls']}|{s['steps']}", []).append(s['device_ms_mean'])
    model = {k: dict(t_ms=statistics.median(v), n=len(v), min=min(v), max=max(v)) for k, v in grp.items()}
    for B in BS:
        for cls in ('const', 'ons'):
            for S in SS: assert f'{B}|{cls}|{S}' in model, ('missing cost cell', B, cls, S)
    return model


def peff(k, n): return (n - 0.5) / n if k == n else k / n


def gpu_pred(t, k, n, B):
    if k == 0: return math.inf, math.inf
    p = peff(k, n); P = 1.0 - (1.0 - p) ** B
    prim = t if P >= 0.995 else t * LN01 / math.log(1.0 - P)
    sec = t * LN01 / (B * math.log(1.0 - p))
    return prim, sec


def main():
    M = cost_model()
    insts = sorted({os.path.basename(p) for p in glob.glob(f'{STUDY}/budget/G*')}, key=lambda g: int(g[1:]))
    labels, cohorts = [], {}
    for ii, g in enumerate(insts):
        for ri, r in enumerate(RULES):
            c = 4 * ii + ri
            for gi, (grid, bd, key) in enumerate([('original', 'budget', 'selected_cfg'), ('extended', 'ext_budget', 'union_selected_cfg')]):
                cands, v80 = [], []
                N = target = BKV = None
                for S in SS:
                    d = json.load(open(f'{STUDY}/{bd}/{g}/{r}_S{S}.json'))
                    if bd == 'budget': N = d['N']
                    fin = d['final']; k, n = fin['k_target'], fin['runs']; target, BKV = fin['target'], fin['BKV']
                    cfg = d[key]; fl = flags(cfg); assert int(cfg['S']) == S
                    for B in BS:
                        t = M[f'{B}|{CLS[r]}|{S}']['t_ms']; prim, sec = gpu_pred(t, k, n, B)
                        cands.append(dict(S=S, B=B, flags=fl, k=k, n=n, t_model_ms=t, pred_primary_ms=prim, pred_secondary_ms=sec))
                    # V80's own secondary-optimal selection with the study's model (round time of 12 engines)
                    tr = fin.get('t_round12_ms')
                    vs = math.inf if (k == 0 or tr is None) else tr * LN01 / (12 * math.log(1.0 - peff(k, n)))
                    v80.append(dict(S=S, flags=fl, k=k, n=n, t_round12_ms=tr, pred_secondary_ms=vs, tts12_ms=fin.get('tts12_ms')))
                for est, field in (('primary', 'pred_primary_ms'), ('secondary', 'pred_secondary_ms')):
                    sel = min(cands, key=lambda x: (x[field], x['S'], x['B']))
                    assert math.isfinite(sel[field]), ('no finite candidate', g, r, grid, est)
                    ck = f"{g}|{sel['flags']}|{sel['B']}"
                    j = 2 * gi + (0 if est == 'primary' else 1)
                    if ck not in cohorts:
                        cohorts[ck] = dict(key=f"{g}_{r.replace('-', '')}_{grid[0]}{est[0]}", instance=g, N=N, target=target, BKV=BKV, rule=r, cls=CLS[r],
                                           S=sel['S'], B=sel['B'], flags=sel['flags'], id_base=110_000_000 + 200_000 * c + 50_000 * j, labels=[])
                    lab = f"{g}|{r}|{grid}|{est}"
                    cohorts[ck]['labels'].append(lab)
                    v80sel = min(v80, key=lambda x: (x['pred_secondary_ms'], x['S']))
                    labels.append(dict(label=lab, instance=g, rule=r, grid=grid, estimator=est, cell=c, cohort=cohorts[ck]['key'], S=sel['S'], B=sel['B'],
                                       flags=sel['flags'], k=sel['k'], n=sel['n'], t_model_ms=sel['t_model_ms'], pred_primary_ms=sel['pred_primary_ms'],
                                       pred_secondary_ms=sel['pred_secondary_ms'], N=N, target=target, BKV=BKV, candidates=cands,
                                       v80_secondary_selection=v80sel if est == 'secondary' else None))
    os.makedirs('results/GS', exist_ok=True)
    json.dump(M, open('results/GS/cost_model.json', 'w'), indent=1)
    json.dump(labels, open('results/GS/plan_labels.json', 'w'), indent=0, default=lambda x: None if x == math.inf else x)
    with open('results/GS/plan_cohorts.tsv', 'w') as f:
        for ck, co in cohorts.items():
            f.write('\t'.join(str(co[x]) for x in ('key', 'instance', 'N', 'target', 'rule', 'cls', 'S', 'B', 'flags', 'id_base')) + '\t' + ','.join(co['labels']) + '\n')
    print(len(labels), 'labels,', len(cohorts), 'unique cohorts')


if __name__ == '__main__':
    main()
