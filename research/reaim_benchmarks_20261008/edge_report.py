"""Grid-edge report of the held-out selections (GPP, TSP): for each method and grid axis, how many of the
(instance, budget) selections sit at the axis minimum or maximum (PROTOCOL.md: "Grid edges will be reported").
Writes edge_report.txt."""
import json
import sys
from collections import defaultdict
from pathlib import Path
HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import grids  # noqa: E402

AX = {'SA': [('T0', 'SA', 'T0'), ('T1', 'SA', 'T1')], 'SCA': [('q', 'SCA', 'q'), ('T0', 'SCA', 'T0'), ('tfin', 'SCA', 'tfin')],
      'TEC': [('jv', 'TEC', 'jv'), ('q', 'shared', 'q'), ('T0', 'shared', 'T0'), ('tfin', 'shared', 'tfin')],
      'Onsager-kT': [('kappa', 'kT', 'kappa'), ('q', 'shared', 'q'), ('T0', 'shared', 'T0'), ('tfin', 'shared', 'tfin')],
      'Onsager-online': [('lam_k2000', 'online', 'lam'), ('q', 'shared', 'q'), ('T0', 'shared', 'T0'), ('tfin', 'shared', 'tfin')],
      'APC-SCA': [('q_reset', 'APC', 'q_reset'), ('r_q', 'APC', 'r_q'), ('q_lim', 'APC', 'q_lim'), ('T0', 'APC', 'T0')],
      'ReAIM ASA': [('kset', None, None), ('T1', None, None)],
      'aSB': [('dt', 'aSB', 'dt'), ('xi', 'aSB', 'xi')], 'bSB': [('dt', 'bSB', 'dt'), ('xi', 'bSB', 'xi')],
      'dSB': [('dt', 'dSB', 'dt'), ('xi', 'dSB', 'xi')]}


def main():
    L = []
    for prob in ('gpp', 'tsp'):
        L.append(f'\n## {prob.upper()} (selections over instances x 5 budgets; counts at axis min / max)\n')
        for meth, axes in AX.items():
            cnt = defaultdict(lambda: [0, 0, 0])
            for f in sorted((HERE / 'results' / prob).glob(f"*/{meth.replace(' ', '_')}_S*.json")):
                d = json.loads(f.read_text()); c = d.get('selected_cfg')
                if not c:
                    continue
                rel = c.get('rel', {})
                for name, gkey, gax in axes:
                    if meth == 'ReAIM ASA':
                        vals = grids.REAIM_KSETS if name == 'kset' else grids.REAIM_T1
                        v = tuple(rel['kset']) if name == 'kset' else rel['T1']
                    elif meth in ('aSB', 'bSB', 'dSB'):
                        vals = grids.G[prob][gkey][gax]; v = c[name]
                    else:
                        vals = grids.G[prob][gkey][gax]; v = rel[name]
                    vals = list(vals) if isinstance(vals, (list, tuple)) else [vals]
                    lo, hi = vals[0] if name != 'T1' or meth != 'ReAIM ASA' else vals[-1], vals[-1] if name != 'T1' or meth != 'ReAIM ASA' else vals[0]
                    if isinstance(v, float):
                        lo_hit = abs(v - min(vals)) < 1e-9 * max(1, abs(v)); hi_hit = abs(v - max(vals)) < 1e-9 * max(1, abs(v))
                    else:
                        lo_hit = v == vals[0]; hi_hit = v == vals[-1]
                    e = cnt[name]; e[0] += 1; e[1] += lo_hit; e[2] += hi_hit and len(vals) > 1
            parts = [f'{k}: {v[1]}/{v[0]} min, {v[2]}/{v[0]} max' for k, v in cnt.items()]
            L.append(f'- {meth}: ' + '; '.join(parts))
    txt = '\n'.join(L)
    (HERE / 'edge_report.txt').write_text(txt + '\n')
    print(txt)


if __name__ == '__main__':
    main()
