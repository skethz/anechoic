"""Addendum 5/5.1: final Table 8b (composition fixed by PROTOCOL_ADDENDUM5_1.md) and the supplementary record.
Every number is read from a final-run file; the source of each row is listed in the output. Writes table8b_a5.md and
table8b_a5.json. Usage: python3 analyze_a5.py"""
import json
import math
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.dont_write_bytecode = True
sys.path.insert(0, str(HERE))
import grids_a2 as GA  # noqa: E402

S_REF = {'mcp': 4000, 'gpp': 4096, 'tsp': 8192}
NINST = {p: len(v) for p, v in GA.INSTANCES.items()}
SEL = json.loads((HERE / 'results_a5' / 'sb_ttt_selection.json').read_text())


def tag(prob, inst):
    return inst if prob == 'tsp' else f'G{inst}'


def a3(m):
    return lambda prob, t, S: HERE / 'results_a3' / 'final' / prob / t / f'{m}_S{S}.json'


def a3bv(m):
    return lambda prob, t, S: HERE / 'results_a3' / 'final_bv' / prob / t / f'{m}_S{S}.json'


def asb09(prob, t, S):
    if prob == 'tsp':
        return HERE / 'results_a4' / 'final' / 'sb_anc' / 'aSB_dt09' / t / f'aSB_S{S}.json'
    return a3('aSB')(prob, t, S)


def asb05(prob, t, S):
    if prob == 'tsp':
        return HERE / 'results_a4' / 'final' / 'sb_anc' / 'aSB_dt05' / t / f'aSB_S{S}.json'
    return HERE / 'results_a4' / 'final' / 'asb_dt05' / prob / t / f'aSB_S{S}.json'


def asbstab(prob, t, S):
    return HERE / 'results_a4' / 'final' / 'asb_stab' / prob / t / f'aSB_S{S}.json'


def reaim_t1(prob, t, S):
    return HERE / 'results_a5' / 'final' / 'reaim_t1' / prob / t / f'ReAIM_ASA_S{S}.json'


def sb_ttt(m):
    def f(prob, t, S):
        dt = SEL[f'{prob}/{t}/{m}']['selected_dt']
        return HERE / 'results_a5' / 'final' / 'sb_ttt' / prob / t / f'{m}_dt{dt:g}_S{S}.json'
    return f


def sb_fixed(m, dt):
    return lambda prob, t, S: HERE / 'results_a5' / 'final' / 'sb_ttt' / prob / t / f'{m}_dt{dt:g}_S{S}.json'


def tec_seq(prob, t, S):
    return HERE / 'results_a4' / 'final' / 'tec_seq' / prob / t / f'TEC_S{S}.json'


MAIN = [  # label, {prob: source function}; composition fixed by PROTOCOL_ADDENDUM5_1.md
    ('SA (Neal defaults)', {p: a3('SA') for p in S_REF}),
    ('STATICA, K2000 setting x σ-transfer [our rule]', {p: a3('SCA') for p in S_REF}),
    ('TEC sequential p-bit, K2000 setting x σ-transfer [our rule]', {p: tec_seq for p in S_REF}),
    ('APC-SCA (GPP: TSP-column settings [our choice])', {p: a3bv('APC-SCA') for p in S_REF}),
    ('ReAIM ASA (Table I k)', {p: reaim_t1 for p in S_REF}),
    ('aSB (Δt 0.9, M 2)', {p: asb09 for p in S_REF}),
    ('bSB (Δt by TTT)', {p: sb_ttt('bSB') for p in S_REF}),
    ('dSB (Δt by TTT)', {p: sb_ttt('dSB') for p in S_REF}),
    ('Onsager-κT (ours)', {p: a3('Onsager-kT') for p in S_REF}),
    ('Onsager-online (ours)', {p: a3('Onsager-online') for p in S_REF}),
]
SUPP = [
    ('aSB Δt 0.5, M 2 [our choice]', {p: asb05 for p in S_REF}),
    ('aSB stability-margin Δt [our rule] (TSP: frozen ancilla)', {p: asbstab for p in S_REF}),
    ('ReAIM k {1,2,6,16} (A3.1) [superseded]', {p: a3bv('ReAIM_ASA') for p in S_REF}),
    ('bSB per-budget Δt by pilot mean quality (A3) [superseded]', {p: a3('bSB') for p in S_REF}),
    ('dSB per-budget Δt by pilot mean quality (A3) [superseded]', {p: a3('dSB') for p in S_REF}),
] + [(f'{m} fixed Δt {dt:g} (A5 finals)', {p: sb_fixed(m, dt) for p in S_REF})
     for m in ('bSB', 'dSB') for dt in (0.25, 0.5, 0.75, 1.0, 1.25)]


def stats(fn, prob):
    insts = GA.INSTANCES[prob]
    q = []; f = []; inv = 0; runs = 0; mcs = {}
    for inst in insts:
        t = tag(prob, inst)
        d = json.loads(fn(prob, t, S_REF[prob]).read_text())['final']
        q.append(d['mean_quality']); f.append(d['p_feasible']); inv += d.get('n_invalid', 0); runs += d['runs']
        mcs[t] = min(float(json.loads(fn(prob, t, S).read_text())['final']['mcs99']) for S in GA.S_LIST[prob])
    solved = [v for v in mcs.values() if math.isfinite(v)]
    return dict(quality=float(np.mean(q)), feasible=float(np.mean(f)), invalid_frac=inv / runs, solved=len(solved),
                n=len(insts), mcs=mcs, geo_own=float(np.exp(np.mean(np.log(solved)))) if solved else None)


def geo(xs):
    return float(np.exp(np.mean(np.log(xs)))) if xs else None


def table(rows):
    T = {}
    for label, spec in rows:
        T[label] = {p: (stats(spec[p], p) if p in spec else None) for p in S_REF}
    return T


def fmt(x, nd=3):
    return '–' if x is None else f'{x:.{nd}f}'


def fs(x):
    return '–' if x is None or not math.isfinite(x) else f'{x:,.0f}'


def cell(e, prob):
    if e is None:
        return ['–'] * 3
    extra = f" (non-finite {100 * e['invalid_frac']:.0f}%)" if e['invalid_frac'] > 0 else ''
    return [fmt(e['quality']) + extra, fmt(e['feasible'], 2), f"{e['solved']}/{e['n']}"]


def main():
    M = table(MAIN); S = table(SUPP)
    insts = [tag('mcp', g) for g in GA.MCP_INST]
    num = [r for r, _ in MAIN]
    sets = {
        'all reported rows': num,
        'all reported rows except aSB': [r for r in num if not r.startswith('aSB')],
        'SA, dSB, Onsager-online (C3)': ['SA (Neal defaults)', 'dSB (Δt by TTT)', 'Onsager-online (ours)'],
    }
    common = {}
    for name, rows in sets.items():
        cs = [t for t in insts if all(math.isfinite(M[r]['mcp']['mcs'][t]) for r in rows)]
        common[name] = dict(instances=cs, geo={r: geo([M[r]['mcp']['mcs'][t] for t in cs]) for r in rows} if cs else {})
    paired = {}
    for ref in ('Onsager-κT (ours)', 'Onsager-online (ours)'):
        for r in num:
            both = [t for t in insts if math.isfinite(M[r]['mcp']['mcs'][t]) and math.isfinite(M[ref]['mcp']['mcs'][t])]
            rat = [M[r]['mcp']['mcs'][t] / M[ref]['mcp']['mcs'][t] for t in both]
            paired[f'{r} / {ref}'] = dict(geo=geo(rat), n=len(both), ref_fewer=sum(1 for x in rat if x > 1))
    c3 = common['SA, dSB, Onsager-online (C3)']
    C3 = dict(n=len(c3['instances']),
              online_over_SA=c3['geo']['Onsager-online (ours)'] / c3['geo']['SA (Neal defaults)'] if c3['geo'] else None,
              online_over_dSB=c3['geo']['Onsager-online (ours)'] / c3['geo']['dSB (Δt by TTT)'] if c3['geo'] else None)
    dtsel = {}
    for key, v in SEL.items():
        prob, _, m = key.split('/')
        dtsel.setdefault(f'{prob}/{m}', {}).setdefault(str(v['selected_dt']), 0)
        dtsel[f'{prob}/{m}'][str(v['selected_dt'])] += 1
    out = dict(main=M, supplement=S, common_sets=common, paired=paired, C3=C3, sb_dt_selected_counts=dtsel)
    (HERE / 'table8b_a5.json').write_text(json.dumps(out, indent=1, default=float) + '\n')

    hdr = ('| Method | Max-Cut quality (S 4000) | Max-Cut feasible | Max-Cut solved/20 | GPP quality (S 4096) | GPP feasible | '
           'GPP solved/9 | TSP quality (S 8192) | TSP valid | TSP solved/6 |')
    L = ['## Table 8b (final, Addendum 5)\n', hdr, '|' + '---|' * 10]
    for r, _ in MAIN:
        L.append(f'| {r} | ' + ' | '.join(sum((cell(M[r][p], p) for p in S_REF), [])) + ' |')
    L += ['\n## Max-Cut steps to the 1% target (MCS99, min over S): geometric means over common solved sets\n']
    for name, c in common.items():
        L.append(f"- {name} ({len(c['instances'])} instances): " +
                 '; '.join(f'{r} {fs(v)}' for r, v in c['geo'].items()))
    L += ['\n## Max-Cut paired step ratios, method / Onsager form, over instances both solve (n; Onsager fewer steps on k)\n',
          '| Method | / Onsager-κT | / Onsager-online |', '|---|---|---|']
    for r in num:
        a = paired[f'{r} / Onsager-κT (ours)']; b = paired[f'{r} / Onsager-online (ours)']
        L.append(f"| {r} | {fmt(a['geo'], 2)}× ({a['n']}, {a['ref_fewer']}) | {fmt(b['geo'], 2)}× ({b['n']}, {b['ref_fewer']}) |")
    L.append(f"\nC3 (instances solved by SA, dSB and Onsager-online: {C3['n']}): Onsager-online/SA = {fmt(C3['online_over_SA'], 3)}, "
             f"Onsager-online/dSB = {fmt(C3['online_over_dSB'], 3)} (ratio of geometric-mean steps; > 1 means more steps for ours)")
    L += ['\n## bSB/dSB: Δt selected by pilot TTT (count of instances)\n']
    for k, v in sorted(dtsel.items()):
        L.append(f'- {k}: ' + ', '.join(f'Δt {a}: {b}' for a, b in sorted(v.items(), key=lambda x: float(x[0]))))
    L += ['\n## Supplementary record (labelled; not Table 8b)\n', hdr, '|' + '---|' * 10]
    for r, _ in SUPP:
        L.append(f'| {r} | ' + ' | '.join(sum((cell(S[r][p], p) for p in S_REF), [])) + ' |')
    (HERE / 'table8b_a5.md').write_text('\n'.join(L) + '\n')
    print('\n'.join(L))


if __name__ == '__main__':
    main()
