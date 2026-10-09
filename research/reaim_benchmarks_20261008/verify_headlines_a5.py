"""Independent check of the Addendum 5/5.1 headline numbers (table8b_a5.json): recomputes, from the per-run arrays of
every final file (quality, values, feasible_runs), the mean quality and feasibility at the ReAIM-matched S, the success
count against the stored target, MCS99 at every S and its minimum, and the instances solved; recomputes the bSB/dSB
time-to-target selection from the pilot files' per-run arrays; and compares everything with table8b_a5.json and
results_a5/sb_ttt_selection.json. Shares no code with analyze_a5.py except the source paths. Usage: python3 verify_headlines_a5.py"""
import json
import math
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.dont_write_bytecode = True
sys.path.insert(0, str(HERE))
import analyze_a5 as A  # noqa: E402  (source paths only)
import grids_a2 as GA  # noqa: E402

DTS = (0.25, 0.5, 0.75, 1.0, 1.25)


def mcs(S, k, n):
    p = k / n
    if p <= 0:
        return math.inf
    if p >= 1:
        return float(S)
    return S * math.log(0.01) / math.log(1 - p)


def run_stats(path, prob):
    d = json.loads(path.read_text())
    q = np.asarray(d['quality'], float); f = np.asarray(d['feasible_runs'], bool); v = np.asarray(d['values'], float)
    t = d['final']['target']
    succ = f & ((v >= t) if prob == 'mcp' else (v <= t))
    return dict(q=float(q.mean()), f=float(f.mean()), k=int(succ.sum()), n=len(q), S=d['S'])


def main():
    T = json.loads((HERE / 'table8b_a5.json').read_text())
    bad = 0; checked = 0
    for label, spec in A.MAIN + A.SUPP:
        for prob in A.S_REF:
            if prob not in spec:
                continue
            fn = spec[prob]; qs = []; fs = []; solved = 0
            for inst in GA.INSTANCES[prob]:
                t = A.tag(prob, inst)
                r = run_stats(fn(prob, t, A.S_REF[prob]), prob)
                qs.append(r['q']); fs.append(r['f'])
                m = min(mcs(S, *[run_stats(fn(prob, t, S), prob)[x] for x in ('k', 'n')]) for S in GA.S_LIST[prob])
                stored = T['main' if (label, spec) in A.MAIN else 'supplement'][label][prob]['mcs'][t]
                checked += 1
                if not ((math.isinf(m) and math.isinf(stored)) or abs(m - stored) <= 1e-6 * max(1.0, m)):
                    bad += 1; print('MCS99 MISMATCH', label, prob, t, m, stored)
                solved += math.isfinite(m)
            e = T['main' if (label, spec) in A.MAIN else 'supplement'][label][prob]
            for name, mine, theirs in (('quality', float(np.mean(qs)), e['quality']), ('feasible', float(np.mean(fs)), e['feasible'])):
                checked += 1
                if abs(mine - theirs) > 1e-9:
                    bad += 1; print('MISMATCH', label, prob, name, mine, theirs)
            checked += 1
            if solved != e['solved']:
                bad += 1; print('SOLVED MISMATCH', label, prob, solved, e['solved'])
    sel = json.loads((HERE / 'results_a5' / 'sb_ttt_selection.json').read_text())
    for prob, insts in GA.INSTANCES.items():
        for inst in insts:
            t = A.tag(prob, inst)
            for m in ('bSB', 'dSB'):
                best = None
                for kk, dt in enumerate(DTS):
                    vals = []
                    for S in GA.S_LIST[prob]:
                        r = run_stats(HERE / 'results_a5' / 'pilot' / 'sb_ttt' / prob / t / f'{m}_dt{dt:g}_S{S}.json', prob)
                        vals.append((mcs(S, r['k'], r['n']), S, r['q']))
                    fin = [x for x in vals if math.isfinite(x[0])]
                    ttt, S_at, q_at = min(fin) if fin else (math.inf, vals[-1][1], vals[-1][2])
                    key = (ttt, -q_at, kk)
                    if best is None or key < best[0]:
                        best = (key, dt)
                checked += 1
                if best[1] != sel[f'{prob}/{t}/{m}']['selected_dt']:
                    bad += 1; print('SELECTION MISMATCH', prob, t, m, best[1], sel[f'{prob}/{t}/{m}']['selected_dt'])
    c3 = T['C3']
    print(f'checked {checked} quantities; mismatches {bad}; C3 n={c3["n"]} online/SA={c3["online_over_SA"]:.4f} online/dSB={c3["online_over_dSB"]:.4f}')
    print('ALL OK' if bad == 0 else 'FAIL')
    return bad == 0


if __name__ == '__main__':
    sys.exit(0 if main() else 1)
