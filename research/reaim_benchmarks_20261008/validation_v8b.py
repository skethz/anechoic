"""Addendum 5 V8b summary (analysis only): ReAIM ASA with Table I k sets against ReAIM Tables IV-VI, with the A3.1
k {1,2,6,16} numbers (validation_a4, V5) beside them. Best of 20 = median over 12 disjoint blocks of 20 runs; TSP: ARPD
over valid tours per block, median over blocks. Writes validation_v8b.json. Usage: python3 validation_v8b.py"""
import json
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
T4 = {'G1': 11624, 'G2': 11615, 'G6': 2168, 'G7': 2006, 'G10': 1991, 'G11': 554, 'G12': 550, 'G13': 576, 'G14': 3046,
      'G19': 900, 'G20': 939}
T5 = {'G1': 7645, 'G2': 7611, 'G3': 7588, 'G4': 7635, 'G5': 7664, 'G14': 1115, 'G15': 1126, 'G16': 1085, 'G17': 1075}
T6 = {'gr17': 16.38, 'gr21': 32.09, 'gr24': 21.84, 'fri26': 43.99, 'bayg29': 31.00, 'bays29': 28.51}
OPT = {'gr17': 2085, 'gr21': 2707, 'gr24': 1272, 'fri26': 937, 'bayg29': 1610, 'bays29': 2020}


def blocks(v, f, how):
    v = np.asarray(v, float); f = np.asarray(f, bool); out = []
    for k in range(12):
        vv = v[20 * k:20 * k + 20][f[20 * k:20 * k + 20]]
        if len(vv):
            out.append(how(vv))
    return (float(np.median(out)) if out else None), 12 - len(out)


def main():
    R = {}
    for name, tab, prob in (('t4', T4, 'mcp'), ('t5', T5, 'gpp'), ('t6', T6, 'tsp')):
        for inst, paper in tab.items():
            row = dict(paper=paper)
            for lab, path in (('table_I_k', HERE / 'validation_a5' / f'reaim_t1_{name}_{inst}.json'),
                              ('k_1_2_6_16', HERE / 'validation_a4' / f'reaim_{name}_{inst}.json')):
                d = json.loads(path.read_text())
                if prob == 'tsp':
                    how = lambda vv, i=inst: float((100 * (vv - OPT[i]) / OPT[i]).mean())  # noqa: E731
                else:
                    how = np.max if prob == 'mcp' else np.min
                val, empty = blocks(d['values'], d['feasible'], how)
                row[lab] = val; row[lab + '_empty_blocks'] = empty; row[lab + '_valid'] = float(np.mean(d['feasible']))
                row[lab + '_dev'] = None if val is None else ((val - paper) if prob == 'tsp' else (val - paper) / paper)
            R[f'{prob}/{inst}'] = row
    (HERE / 'validation_v8b.json').write_text(json.dumps(R, indent=1) + '\n')
    for k, r in R.items():
        if k.startswith('tsp'):
            print(f"{k:12s} paper ARPD {r['paper']:6.2f} | Table I k: {r['table_I_k']} (valid {r['table_I_k_valid']:.2f}) | "
                  f"k{{1,2,6,16}}: {r['k_1_2_6_16']:.2f} (valid {r['k_1_2_6_16_valid']:.2f})")
        else:
            print(f"{k:12s} paper {r['paper']:6d} | Table I k: {r['table_I_k']:8.1f} ({100 * r['table_I_k_dev']:+.2f}%) | "
                  f"k{{1,2,6,16}}: {r['k_1_2_6_16']:8.1f} ({100 * r['k_1_2_6_16_dev']:+.2f}%)")


if __name__ == '__main__':
    main()
