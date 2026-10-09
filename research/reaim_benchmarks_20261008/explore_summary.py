"""Tables of the exploratory (post hoc) runs: GPP penalty sensitivity, ReAIM F = min, and SA/bSB K-bit sensitivity.
Writes explore/explore_summary.txt."""
import json
from pathlib import Path
import numpy as np
HERE = Path(__file__).resolve().parent
E = HERE / 'explore'
L = []
p = E / 'penalty'
if p.exists():
    L.append('## GPP penalty sensitivity (P = 4A; protocol P = 4), S = 4096: mean quality (feasible) [P(target)]\n')
    recs = [json.loads(f.read_text()) for f in sorted(p.glob('*/*.json'))]
    meths = ['SA', 'bSB', 'SCA', 'TEC', 'APC-SCA', 'Onsager-kT', 'Onsager-online', 'ReAIM ASA']
    L.append('| Instance, P | ' + ' | '.join(meths) + ' |'); L.append('|---|' + '---|' * len(meths))
    prot = {}
    for g in (1, 14):
        for m in meths:
            f = HERE / 'results' / 'gpp' / f'G{g}' / f"{m.replace(' ', '_')}_S4096.json"
            if f.exists():
                d = json.loads(f.read_text())['final']; prot[(f'G{g}', 4, m)] = (d['mean_quality'], d['p_feasible'], d['p_target'])
    for r in recs:
        prot[(r['instance'], r['P'], r['method'])] = (r['mean_quality'], r['p_feasible'], r['p_target'])
    for g in ('G1', 'G14'):
        for P in (1, 2, 4, 8):
            cells = [('%.3f (%.2f) [%.2f]' % prot[(g, P, m)]) if (g, P, m) in prot else '–' for m in meths]
            L.append(f'| {g}, P={P}{" (protocol)" if P == 4 else ""} | ' + ' | '.join(cells) + ' |')
p = E / 'reaim_F'
if p.exists():
    L.append('\n## ReAIM ASA with F = min (post hoc) vs the protocol F = max, at ReAIM budgets\n')
    L.append('| Instance | F=max quality (feas) [P(target)] | F=min quality (feas) [P(target)] | Step-4 pick | F=min ARPD (feasible) |')
    L.append('|---|---|---|---|---|')
    for f in sorted(p.glob('*/*.json')):
        r = json.loads(f.read_text()); pm = r['protocol_F_max']
        L.append(f"| {r['instance']} | {pm['final_mean_quality']:.3f} ({pm['final_p_feasible']:.2f}) [{pm['final_p_target']:.2f}] | "
                 f"{r['mean_quality']:.3f} ({r['p_feasible']:.2f}) [{r['p_target']:.2f}] | {r['step4_choice']} | "
                 f"{r.get('arpd_feasible', float('nan')):.1f} |")
p = E / 'precision_others'
if p.exists():
    L.append('\n## SA and bSB on K-bit quantized GPP/TSP (float model, post hoc): mean quality (feasible)\n')
    L.append('| Instance | Method | full | K=2 | K=3 | K=4 | K=6 | K=8 |'); L.append('|---|---|---|---|---|---|---|---|')
    for f in sorted(p.glob('*.json')):
        r = json.loads(f.read_text()); b = r['by_K']
        L.append(f"| {r['instance']} | {r['method']} | " + ' | '.join('%.3f (%.2f)' % (b[k]['mean_quality'], b[k]['p_feasible'])
                                                            for k in ('None', '2', '3', '4', '6', '8')) + ' |')
txt = '\n'.join(L)
(E / 'explore_summary.txt').write_text(txt + '\n')
print(txt)
