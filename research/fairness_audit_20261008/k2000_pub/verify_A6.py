"""Independent check of PROTOCOL_PUB_A6 outputs (shares no code with sb_ttt_A6.py, pub_methods.py or run_pub.py):
  1. the selected Delta t, recomputed from the raw pilot cuts (TTT = min over S of S ln0.01/ln(1-p), p = P(cut >= 33000);
     ties -> higher pilot mean cut at that S, then smaller Delta t);
  2. results_pub_A6/S*.json: bSB/dSB entries equal the raw finals at the selected Delta t, and their mean_cut, p33000 and
     mcs99 recomputed from the stored per-run cuts; every other entry identical to results_pub_A2/S*.json;
  3. the pre-declared kernel check: A6 finals at Delta t = 1, S = 1000 against the COP study's V7 runs (solvers.sb)
     within 3 combined SE, and against Goto 2021 Fig. 2A read-offs within 40 + 3 SE.
Usage: python3 verify_A6.py"""
import json
import math
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
OUT = HERE / 'results_pub_A6'
A2 = HERE / 'results_pub_A2'
V7 = HERE.parents[1] / 'reaim_benchmarks_20261008' / 'validation_a4'
S_LIST = (250, 500, 1000, 2000, 4000)
DTS = (0.25, 0.5, 0.75, 1.0, 1.25)
GOTO_1000 = {'bSB': 33215, 'dSB': 33110}


def mcs(S, p):
    if p <= 0:
        return math.inf
    if p >= 1:
        return float(S)
    return S * math.log(0.01) / math.log(1 - p)


def raw(name, S, k, kind):
    return json.loads((OUT / 'raw' / f'{name}_S{S}_k{k}_{kind}.json').read_text())


def main():
    bad = 0
    sel = json.loads((OUT / 'selection_A6.json').read_text())
    chosen = {}
    for name in ('bSB', 'dSB'):
        keys = []
        for k, dt in enumerate(DTS):
            vals = []
            for S in S_LIST:
                c = np.array(raw(name, S, k, 'pilot')['cuts'], float)
                assert len(c) == 256
                vals.append((mcs(S, float((c >= 33000).mean())), S, float(c.mean())))
            t, S_at, mc = min(vals, key=lambda v: (v[0], v[1]))
            keys.append(((t, -mc, k), dt))
        dt = min(keys)[1]; chosen[name] = DTS.index(dt)
        ok = dt == sel[name]['selected_dt']; bad += not ok
        print(f'{name}: independent selection dt {dt} vs selection_A6 {sel[name]["selected_dt"]} -> {"OK" if ok else "MISMATCH"}')
    for S in S_LIST:
        a2 = json.loads((A2 / f'S{S}.json').read_text()); a6 = json.loads((OUT / f'S{S}.json').read_text())
        if list(a2) != list(a6):
            bad += 1; print('key order differs', S)
        for n in a2:
            if n in ('bSB', 'dSB'):
                f = a6[n]['final']; r = raw(n, S, chosen[n], 'final')
                same = json.dumps(f) == json.dumps(r)
                c = np.array(f['cuts'], float); p = float((c >= 33000).mean())
                rec = (abs(c.mean() - f['mean_cut']) < 1e-9 and abs(p - f['p33000']) < 1e-12 and
                       (math.isinf(mcs(S, p)) and math.isinf(f['mcs99']) or abs(mcs(S, p) - f['mcs99']) < 1e-9) and
                       len(f['H_mean']) == S + 1 and f['job']['cfg']['dt'] == DTS[chosen[n]] and f['runs'] == 256)
                if not (same and rec):
                    bad += 1; print('bSB/dSB entry check failed', S, n, same, rec)
            elif json.dumps(a2[n]) != json.dumps(a6[n]):
                bad += 1; print('entry changed', S, n)
    for name in ('bSB', 'dSB'):
        r = raw(name, 1000, 3, 'final'); c = np.array(r['cuts'], float)
        v = []
        for k in range(4):
            v += json.loads((V7 / f'sb21_{name}_S1000_c{k}.json').read_text())['cuts']
        v = np.array(v, float)
        se = math.sqrt(c.var(ddof=1) / len(c) + v.var(ddof=1) / len(v)); se6 = c.std(ddof=1) / math.sqrt(len(c))
        ok1 = abs(c.mean() - v.mean()) <= 3 * se; ok2 = abs(c.mean() - GOTO_1000[name]) <= 40 + 3 * se6
        bad += (not ok1) + (not ok2)
        print(f'{name} dt 1 S 1000: A6 {c.mean():.1f} (SE {se6:.2f}) vs V7 {v.mean():.1f} -> diff {c.mean() - v.mean():+.1f}, '
              f'3 SE {3 * se:.1f}: {"OK" if ok1 else "FAIL"}; vs Goto Fig. 2A ~{GOTO_1000[name]}: {"OK" if ok2 else "FAIL"}')
    print('ALL OK' if bad == 0 else f'FAIL ({bad})')


if __name__ == '__main__':
    main()
