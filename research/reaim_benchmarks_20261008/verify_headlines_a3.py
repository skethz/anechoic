"""Independent check of the A3 headline numbers (Table 8b) before reporting. Writes verify_headlines_a3.json.
 1. From the raw per-run arrays of every final file, recompute quality, feasibility, P(target) and MCS99 with separate
    code (G-set best-known values from bkv.json, GPP references from gpp_reference.json, TSP optima from the TSPLIB
    list) and compare with the stored statistics.
 2. Re-run a sample of finals (one per method and problem, at the smallest budget) with the stored configuration and the
    recomputed seed, decode the spins with independent decoders (G-set edge file / TSPLIB matrix), and compare the
    per-run objective values with the stored ones.
 3. Recompute the Table 8b headline numbers (mean quality and feasibility at the ReAIM budget, common-set Max-Cut step
    geo-means and the SA/dSB vs Onsager-online ratios) and compare with results_summary_a3.json."""
import json
import math
import os
import sys
from pathlib import Path

sys.dont_write_bytecode = True
HERE = Path(__file__).resolve().parent
os.environ.setdefault('NUMBA_CACHE_DIR', str(HERE / '.numba_cache'))
sys.path.insert(0, str(HERE))
import numpy as np  # noqa: E402

GSET = HERE.parent / 'gset_20261007'
BKV = json.loads((GSET / 'bkv.json').read_text())['instances']
REF = json.loads((HERE / 'data' / 'gpp_reference.json').read_text())['instances']
OPT = {'gr17': 2085, 'gr21': 2707, 'gr24': 1272, 'fri26': 937, 'bayg29': 1610, 'bays29': 2020}
METHODS = ('SA', 'SCA', 'TEC', 'APC-SCA', 'ReAIM ASA', 'aSB', 'bSB', 'dSB', 'Onsager-kT', 'Onsager-online')
S_LIST = {'mcp': (250, 500, 1000, 2000, 4000), 'gpp': (256, 512, 1024, 2048, 4096), 'tsp': (512, 1024, 2048, 4096, 8192)}
INST = {'mcp': [f'G{g}' for g in range(1, 21)], 'gpp': [f'G{g}' for g in (1, 2, 3, 4, 5, 14, 15, 16, 17)],
        'tsp': ['gr17', 'gr21', 'gr24', 'fri26', 'bayg29', 'bays29']}
S_REAIM = {'mcp': 4000, 'gpp': 4096, 'tsp': 8192}


def load(prob, t, m, S):
    return json.loads((HERE / 'results_a3' / 'final' / prob / t / f"{m.replace(' ', '_')}_S{S}.json").read_text())


def per_run(prob, t, d):
    v = np.array(d['values'], float); f = np.array(d['feasible_runs'], bool)
    if prob == 'mcp':
        b = BKV[t]; q = v / b['BKV']; succ = v >= b['target']
    elif prob == 'gpp':
        q = np.where(f, REF[t]['R'] / np.maximum(v, 1), 0.0); succ = f & (v <= REF[t]['target'])
    else:
        q = np.where(f, OPT[t] / np.maximum(v, 1), 0.0); succ = f & (v <= math.floor(1.01 * OPT[t]))
    if d.get('finite') is not None:
        ok = np.array(d['finite'], bool); q = np.where(ok, q, 0.0); succ &= ok; f &= ok
    return q, f, succ


def mcs99(S, p):
    return math.inf if p <= 0 else (float(S) if p >= 1 else S * math.log(0.01) / math.log1p(-p))


def gset_edges(t):
    lines = (GSET / 'data' / t).read_text().split('\n')
    a = np.array([ln.split() for ln in lines[1:] if ln.strip()], dtype=np.int64)
    return a[:, 0] - 1, a[:, 1] - 1, a[:, 2]


def decode(prob, t, s):
    s = np.asarray(s)
    if prob in ('mcp', 'gpp'):
        i, j, w = gset_edges(t)
        cut = ((s[:, i] != s[:, j]) * w[None, :]).sum(1)
        if prob == 'mcp':
            return cut.astype(float)
        bal = s.sum(1) == 0
        return np.where(bal, cut, -1).astype(float)
    import tsplib
    W = tsplib.load(t)['W']; n = len(W)
    out = []
    for r in s:
        x = (r.reshape(n, n) > 0)
        if not (np.all(x.sum(0) == 1) and np.all(x.sum(1) == 1)):
            out.append(-1.0); continue
        tour = [int(np.argmax(x[:, j])) for j in range(n)]
        out.append(float(tsplib.tour_length(W, tour)))
    return np.array(out)


def main():
    res = dict(stat_mismatches=[], rerun=[], headline={}); ok = True
    for prob in ('mcp', 'gpp', 'tsp'):
        for t in INST[prob]:
            for m in METHODS:
                for S in S_LIST[prob]:
                    d = load(prob, t, m, S); fi = d['final']
                    q, f, succ = per_run(prob, t, d)
                    p = succ.mean()
                    checks = dict(q=abs(q.mean() - fi['mean_quality']) < 1e-12, f=abs(f.mean() - fi['p_feasible']) < 1e-12,
                                  p=abs(p - fi['p_target']) < 1e-12, mcs=(mcs99(S, p) == fi['mcs99']) or
                                  (math.isfinite(mcs99(S, p)) and abs(mcs99(S, p) - fi['mcs99']) < 1e-6))
                    if not all(checks.values()):
                        ok = False; res['stat_mismatches'].append(dict(prob=prob, inst=t, method=m, S=S, checks=checks))
    print('stat mismatches:', len(res['stat_mismatches']), flush=True)
    # 2. re-run sample
    import run_a2 as RA
    import solvers_a3 as SV3
    pcode = {'mcp': 0, 'gpp': 1, 'tsp': 2}
    sample = {'mcp': 'G14', 'gpp': 'G1', 'tsp': 'gr17'}
    for prob, t in sample.items():
        inst = t if prob == 'tsp' else int(t[1:])
        P = RA.problem(prob, inst)
        S = S_LIST[prob][0]
        for m in METHODS:
            d = load(prob, t, m, S)
            mi = METHODS.index(m); si = 0
            s, _, finite = SV3.run_method(P, d['cfg'], 256, np.random.SeedSequence([20261008, 41, pcode[prob], RA.icode(prob, inst), mi, si]))
            v = decode(prob, t, s)
            stored = np.array(d['values'], float)
            if finite is not None:
                okf = np.array(finite, bool); v = np.where(okf, v, stored)  # invalid runs: compare only the finite ones
            same = bool(np.array_equal(v, stored))
            ok &= same
            res['rerun'].append(dict(prob=prob, inst=t, method=m, S=S, identical=same))
            print('rerun', prob, t, m, S, same, flush=True)
    # 3. headline numbers
    summ = json.loads((HERE / 'results_summary_a3.json').read_text())['summary']
    for prob in ('mcp', 'gpp', 'tsp'):
        S = S_REAIM[prob]
        for m in METHODS:
            qs = []; fs = []
            for t in INST[prob]:
                q, f, _ = per_run(prob, t, load(prob, t, m, S)); qs.append(q.mean()); fs.append(f.mean())
            mine = (float(np.mean(qs)), float(np.mean(fs)))
            theirs = (summ[prob]['rows'][m]['q3'], summ[prob]['rows'][m]['f3'])
            same = abs(mine[0] - theirs[0]) < 1e-12 and abs(mine[1] - theirs[1]) < 1e-12
            ok &= same
            res['headline'][f'{prob}/{m}'] = dict(quality=mine[0], feasible=mine[1], matches_analysis=bool(same))
    # Max-Cut steps on the common solved set
    mins = {m: {} for m in METHODS}
    for m in METHODS:
        for t in INST['mcp']:
            best = math.inf
            for S in S_LIST['mcp']:
                _, _, succ = per_run('mcp', t, load('mcp', t, m, S)); best = min(best, mcs99(S, succ.mean()))
            mins[m][t] = best
    common = [t for t in INST['mcp'] if all(math.isfinite(mins[m][t]) for m in METHODS)]
    g = {m: math.exp(np.mean([math.log(mins[m][t]) for t in common])) if common else None for m in METHODS}
    res['headline']['mcp_common_set'] = common
    res['headline']['mcp_steps_common_geo'] = g
    res['headline']['mcp_solved'] = {m: sum(1 for t in INST['mcp'] if math.isfinite(mins[m][t])) for m in METHODS}
    if common:
        res['headline']['online_over_SA'] = g['Onsager-online'] / g['SA']
        res['headline']['online_over_dSB'] = g['Onsager-online'] / g['dSB']
        mine_g = [summ['mcp']['rows'][m]['geo_common3'] for m in METHODS]
        same = all(abs(a - b) < 1e-6 * max(1, abs(b)) for a, b in zip([g[m] for m in METHODS], mine_g))
        ok &= same; res['headline']['common_geo_matches_analysis'] = bool(same)
    res['all_ok'] = bool(ok)
    (HERE / 'verify_headlines_a3.json').write_text(json.dumps(res, indent=1, default=str) + '\n')
    print(json.dumps(res['headline'], indent=1, default=str))
    print('ALL OK' if ok else 'MISMATCH')


if __name__ == '__main__':
    main()
