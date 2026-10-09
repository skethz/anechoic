#!/usr/bin/env python3
"""Phase-G run plan (PROTOCOL.md section 5) from the software study's selected configurations.
Input : data/gset_selected_configs.json (copy of research/gset_20261007/selected_configs.json, sha256 in the protocol)
Output: results/G/plan.json and results/G/plan.tsv (one line per cell: tag, instance, N, target, kind, rule, class, S, flags,
        trial-id base). Deterministic; no randomness.
Mapping of the study's rules onto the engine's V80 host table flags (same table code, tables.hpp):
  plain   -> --lambda 0                    TEC        -> --tec-jv jv
  tecT    -> --tecT kappa --ramp           onsager    -> --lambda lam --ramp
  common  -> --q q --t0 T0 --t1 tfin --steps S   (absolute units of the G-set couplings, J = -w)
APC-SCA is not an engine mode and is omitted."""
import json, os, sys

RULES = ['SCA', 'TEC', 'Onsager-kT', 'Onsager-online']
VERSION = 'extended_amendment1'


def flags(cfg):
    f = cfg['family']
    common = f"--q {cfg['q']!r} --t0 {cfg['T0']!r} --t1 {cfg['tfin']!r} --steps {int(cfg['S'])}"
    if f == 'plain':
        return '--lambda 0 ' + common, 'const'
    if f == 'tec':
        return f"--tec-jv {cfg['jv']!r} " + common, 'const'
    if f == 'tecT':
        assert cfg.get('ramp') is True
        return f"--tecT {cfg['kappa']!r} --ramp " + common, 'const'
    if f == 'onsager':
        assert cfg.get('ramp') is True
        return f"--lambda {cfg['lam']!r} --ramp " + common, 'ons'
    raise ValueError(f)


d = json.load(open('data/gset_selected_configs.json'))
inst = sorted(d['instances'], key=lambda g: int(g[1:]))
cells = []
for ii, g in enumerate(inst):
    v = d['instances'][g]
    for ri, r in enumerate(RULES):
        c = 4 * ii + ri
        sel = v['rules'][r][VERSION]
        for kind in ['E12', 'E1']:
            t = sel['tts_' + kind]
            fl, cls = flags(t['cfg'])
            assert int(t['cfg']['S']) == int(t['S'])
            cells.append(dict(cell=c, tag=f"{g}_{r.replace('-', '')}_{kind}", instance=g, N=v['N'], target=v['target'], BKV=v['BKV'], kind=kind,
                              rule=r, family=t['cfg']['family'], cls=cls, S=int(t['S']), flags=fl, id_base=40_000_000 + 100_000 * c,
                              study_tts_ms=t['tts_ms'], study_p=t['p']))
os.makedirs('results/G', exist_ok=True)
json.dump(cells, open('results/G/plan.json', 'w'), indent=0)
with open('results/G/plan.tsv', 'w') as f:
    for x in cells:
        f.write('\t'.join(str(x[k]) for k in ['tag', 'instance', 'N', 'target', 'kind', 'rule', 'cls', 'S', 'flags', 'id_base']) + '\n')
print(len(cells), 'cells,', len(inst), 'instances')
