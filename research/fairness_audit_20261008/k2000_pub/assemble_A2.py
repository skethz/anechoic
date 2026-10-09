"""PROTOCOL_PUB_A2.md: Figure 2 data with the COP study's readings, from the existing PROTOCOL_PUB runs (no new samples):
STATICA T_init 40 at every budget; APC-SCA and ReAIM with their papers' output (best visited / x_best); bSB and dSB with the
COP reading of Goto's Delta t rule (per budget, highest pilot mean cut, ties to the smaller Delta t). TEC and aSB: from
the COP study's VALIDATED_BASELINES.md when available (argument), else the PROTOCOL_PUB TEC and the A1 aSB (provisional).
Also prints the claim checks. Usage: python3 assemble_A2.py [cop_tec_asb.json]"""
import json
import math
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
RES = HERE / 'results_pub'
OUT = HERE / 'results_pub_A2'
ORIG = HERE.parents[1] / 'algorithm_compare_20261007/results'
S_LIST = (250, 500, 1000, 2000, 4000)
SB = ('aSB', 'bSB', 'dSB')
OURS = ('TEC-T (ours)', 'Onsager SCA (ours)')


def raw(name, S, k, kind='final'):
    tag = name.replace(' ', '_').replace('(', '').replace(')', '')
    return json.loads((RES / 'raw' / f'{tag}_S{S}_k{k}_{kind}.json').read_text())


def mcs99(S, p):
    return math.inf if p <= 0 else (S if p >= 1 else S * math.log(0.01) / math.log1p(-p))


def main():
    plot = {S: json.loads((RES / f'S{S}.json').read_text()) for S in S_LIST}
    a1 = json.loads((RES / 'asb_perrun_A1.json').read_text())
    xb = json.loads((RES / 'reaim_xbest_A2.json').read_text()) if (RES / 'reaim_xbest_A2.json').exists() else None
    cop = json.loads(Path(sys.argv[1]).read_text()) if len(sys.argv) > 1 else None
    out = {S: {} for S in S_LIST}; notes = {}
    for S in S_LIST:
        P = plot[S]
        out[S]['SA'] = dict(final=P['SA (Neal)']['final'], setting='Neal defaults')
        f = raw('SCA (STATICA)', S, 1); out[S]['SCA'] = dict(final=f, setting='STATICA q 4, T 40 -> 5')
        # APC: paper output = argmin over s including s = 1 (best_visited)
        f = dict(P['APC-SCA (published)']['final']); bv = f['best_visited']
        f.update(p33000=bv['p33000'], p33000_wilson95=bv['p33000_wilson95'], k33000=bv['k33000'], mcs99=bv['mcs99'],
                 mean_cut=bv['mean_cut'], output='argmin_s H (paper)')
        out[S]['APC-SCA'] = dict(final=f, setting='APC Alg. 2, r_q 0.45, q_limit 0, T 10 -> 0.1')
        f = dict(P['ReAIM ASA']['final'])
        if xb:
            x = xb[str(S)]; f.update(p33000=x['p33000'], p33000_wilson95=x['p33000_wilson95'], k33000=x['k33000'],
                                      mcs99=x['mcs99'], mean_cut=x['mean_cut'], cuts=x['cuts_xbest'], output='x_best (paper)')
        out[S]['ReAIM ASA'] = dict(final=f, setting='ReAIM T 1 -> 0.1, F max, k {128..1024}')
        # bSB / dSB: COP reading (per budget, highest pilot mean cut, ties -> smaller dt)
        for name in ('bSB', 'dSB'):
            pil = [np.mean(raw(name, S, k, 'pilot')['cuts']) for k in range(5)]
            k = max(range(5), key=lambda i: (pil[i], -i))
            out[S][name] = dict(final=raw(name, S, k), setting=f'Goto 2021, dt {[0.25, 0.5, 0.75, 1.0, 1.25][k]} (COP reading)')
            notes[(name, S)] = k
        if cop and str(S) in cop.get('aSB', {}):
            out[S]['aSB'] = dict(final=cop['aSB'][str(S)], setting='COP VALIDATED_BASELINES')
        else:
            x = a1[str(S)]; f = dict(P['aSB']['final'])
            f.update(p33000=x['p33000'], p33000_wilson95=x['p33000_wilson95'], k33000=x['k33000'], mcs99=x['mcs99'],
                     mean_cut=x['mean_cut_nondiverged'] if x['mean_cut_nondiverged'] is not None else float('-inf'),
                     H_mean=x['H_mean_nondiverged'] or f['H_mean'], diverged=x['diverged'])
            out[S]['aSB'] = dict(final=f, setting='Goto 2019 dt 0.9, M 2 (A1 per run; provisional)')
        if cop and str(S) in cop.get('TEC', {}):
            out[S]['TEC'] = dict(final=cop['TEC'][str(S)], setting='COP VALIDATED_BASELINES')
        else:
            out[S]['TEC'] = dict(final=P['TEC (published)']['final'], setting='TEC Glauber SA, J_v 30, T 100 -> 0.1 (provisional)')
        for n in OURS:
            out[S][n] = P[n]
    OUT.mkdir(exist_ok=True)
    for S in S_LIST:
        (OUT / f'S{S}.json').write_text(json.dumps(out[S]) + '\n')
    names = list(out[S_LIST[0]].keys())
    L = ['A2 (COP-aligned readings): mean cut / p33000 / MCS99 per budget']
    for S in S_LIST:
        L.append(f'--- S = {S}')
        for n in sorted(names, key=lambda n: -out[S][n]['final']['mean_cut']):
            f = out[S][n]['final']
            L.append(f"{n:20s} mean {f['mean_cut']:9.1f}  p {f['p33000']:.3f}  MCS99 {f['mcs99']:9.0f}   [{out[S][n].get('setting', 'frozen')}]")
    L.append('Claim (a) lowest mean energy among discrete methods (all but aSB, bSB, dSB):')
    for S in S_LIST:
        disc = [n for n in names if n not in SB and n not in OURS]
        bb = max(disc, key=lambda n: out[S][n]['final']['mean_cut'])
        ok = all(out[S][o]['final']['mean_cut'] > out[S][bb]['final']['mean_cut'] for o in OURS)
        L.append(f"  S={S}: best discrete baseline {bb} {out[S][bb]['final']['mean_cut']:.1f} -> both forms lower energy: {ok}")
    L.append('Steps to solution (min over S):')
    mins = {n: min((out[S][n]['final']['mcs99'], S) for S in S_LIST) for n in names}
    for n in sorted(names, key=lambda n: mins[n][0]):
        L.append(f'  {n:20s} {mins[n][0]:9.0f} at S = {mins[n][1]}')
    L.append('S = 500 success: ' + ', '.join(f"{n} {out[500][n]['final']['p33000']:.3f}" for n in names))
    L.append(f'bSB/dSB per-budget dt index (COP reading): {notes}')
    t = '\n'.join(L)
    print(t)
    (OUT / 'claims_A2.txt').write_text(t + '\n')


if __name__ == '__main__':
    main()
