"""Held-out analysis for PROTOCOL_ABL.md: device TTS per family and E, hypotheses A1-A4."""
import json
import math
from pathlib import Path

import abl

ROOT = Path(__file__).resolve().parent
d = json.loads((ROOT / 'holdout.json').read_text())
sel, hold = d['selection'], d['holdout']
KEYS = ('family', 'q', 'T0', 'S', 'lam', 'ramp', 'jv', 'kappa', 'j_lo', 'j_hi', 'q_reset', 'r_q', 'q_lim')


def key(c):
    return tuple((k, c.get(k)) for k in KEYS)


by = {key(r): r for r in hold if r['family'] != 'onsager_table'}
tables = {key(dict(r, family='onsager')): r for r in hold if r['family'] == 'onsager_table'}
fams = ['plain', 'tec', 'tecT', 'tecS', 'apc', 'onsager']
res = {}
print(f"{'family':8s} " + " ".join(f"{'E=' + str(E):>26s}" for E in abl.E_LIST))
for fam in fams:
    row = []
    for E in abl.E_LIST:
        c = sel[f'{fam}_E{E}']; r = by[key(c)]
        tts = abl.dev_tts(r['t_ms'], r['p'], E)
        res[(fam, E)] = dict(cfg=c, p=r['p'], t_ms=r['t_ms'], tts=tts)
        row.append(f"p={r['p']:.3f} TTS={tts:7.3f} ms")
    print(f"{fam:8s} " + " ".join(f"{x:>26s}" for x in row))
print()
verdict = {}
for name, fam in (('A1', 'tecT'), ('A2', 'tecS'), ('A3', 'apc')):
    ok = True; parts = []
    for E in abl.E_LIST:
        ratio = res[(fam, E)]['tts'] / res[('onsager', E)]['tts']
        parts.append(f"E={E}: {ratio:.2f}x"); ok &= ratio >= 1.2
    verdict[name] = ok
    print(f"{name} Onsager vs {fam}: " + ", ".join(parts) + f" -> {'PASS' if ok else 'FAIL'} (need >= 1.2x at every E)")
print("A4 online vs fixed table (same configs, paired):")
for k, r in tables.items():
    on = by[k]
    print(f"   {dict(k)['T0']=}, S={dict(k)['S']}, lam={dict(k)['lam']}: online p={on['p']:.3f} flips={on['mean_flips']:.0f} | "
          f"table p={r['p']:.3f} flips={r['mean_flips']:.0f}")
(ROOT / 'summary.json').write_text(json.dumps(dict(verdict=verdict, results={f'{f}_E{E}': v for (f, E), v in res.items()}), indent=1) + '\n')
