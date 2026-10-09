"""PROTOCOL_PUB_A6.md: K2000 bSB and dSB with ONE Delta t per method chosen by time-to-target (Goto et al. 2021), the
rule of the COP study's Table 8b (research/reaim_benchmarks_20261008/PROTOCOL_ADDENDUM5.md, A5.2), so that Figure 2 and
Table 8b agree. Kernel, budgets, run counts and summary format are those of PROTOCOL_PUB (methods.sb_traj via
pub_methods.run_vectorized; pub_methods.summarize_final); only the seeds (fresh) and the Delta t rule differ from A2.

Jobs (256 runs each; S in {250, 500, 1000, 2000, 4000}; Delta t index k over {0.25, 0.5, 0.75, 1, 1.25}):
  pilot  SeedSequence([20261008, 701, code, S index, k])  -> results_pub_A6/raw/<m>_S<S>_k<k>_pilot.json
  final  SeedSequence([20261008, 702, code, S index, k])  -> results_pub_A6/raw/<m>_S<S>_k<k>_final.json
Selection ('select'): run_pub.sb_select on the pilots (PROTOCOL_PUB's frozen rule: minimise over Delta t the minimum
over S of S ln(0.01)/ln(1 - p), p = P(cut >= 33,000); ties -> higher pilot mean cut at that S, then smaller Delta t).
Sensitivity (reported, not used): target 33,004 = ceil(0.99 x 33,337) (Goto's "99% of the best known"), and Goto's
TTT convention that caps the step count at T_com (= S) when P_S > 0.99.
Assembly ('assemble'): results_pub_A6/S<S>.json = results_pub_A2/S<S>.json with the bSB and dSB entries replaced by
the A6 finals at the selected Delta t (every other entry copied unchanged); claims_A6.txt; compare_A6_vs_A2.json.
Usage: python3 sb_ttt_A6.py check | run [workers] | select | assemble"""
import hashlib
import json
import math
import os
import sys
import time
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path

import numpy as np

sys.dont_write_bytecode = True
HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import pub_methods as P  # noqa: E402
import run_pub as RP  # noqa: E402  (sb_select only; no runner code is called)
from common import M  # noqa: E402

OUT = HERE / 'results_pub_A6'
RAW = OUT / 'raw'
A2 = HERE / 'results_pub_A2'
S_LIST = (250, 500, 1000, 2000, 4000)
B = 256
SEED = 20261008
METHODS = ('bSB', 'dSB')
assert RP.S_LIST == S_LIST and tuple(P.SB_DTS) == (0.25, 0.5, 0.75, 1.0, 1.25)


def jobs():
    out = []
    for si, S in enumerate(S_LIST):
        for name in METHODS:
            for k in range(len(P.SB_DTS)):
                for kind, fam in (('pilot', 701), ('final', 702)):
                    out.append(dict(name=name, S=S, si=si, k=k, kind=kind, cfg=P.cfg_of(name, S, k),
                                    seed=[SEED, fam, P.CODE[name], si, k], chunk=None))
    return out


def job_path(j):
    return RAW / f"{j['name']}_S{j['S']}_k{j['k']}_{j['kind']}.json"


_J = None


def _init():
    global _J
    np.seterr(all='ignore')
    _J = M.m.load()


def _run(j):
    J, sumw = _J; t0 = time.time()
    E, cuts, _ = P.run_vectorized(J, sumw, j['cfg'], B, np.random.SeedSequence(j['seed']))
    if j['kind'] == 'pilot':
        res = dict(cuts=np.asarray(cuts).tolist(), mean_flips=None)
    else:
        res = P.summarize_final(E, cuts, sumw, j['S'])
    res.update(job=dict(j), sec=time.time() - t0)
    path = job_path(j); tmp = path.with_suffix('.tmp')
    tmp.write_text(json.dumps(res) + '\n'); tmp.rename(path)
    return path.name, res['sec']


def run(workers):
    RAW.mkdir(parents=True, exist_ok=True)
    todo = [j for j in jobs() if not job_path(j).exists()]
    todo.sort(key=lambda j: -j['S'])
    print(f'A6: {len(todo)} jobs, {workers} workers, start {time.strftime("%Y-%m-%d %H:%M:%S")}', flush=True)
    t0 = time.time()
    with ProcessPoolExecutor(max_workers=workers, initializer=_init) as ex:
        futs = [ex.submit(_run, j) for j in todo]
        for n, f in enumerate(as_completed(futs), 1):
            name, sec = f.result()
            print(f'[{n}/{len(todo)}] {time.time() - t0:6.0f}s {name} ({sec:.0f}s)', flush=True)
    print(f'A6 complete {time.strftime("%Y-%m-%d %H:%M:%S")} ({time.time() - t0:.0f}s)', flush=True)


def load(name, S, k, kind):
    return json.loads((RAW / f'{name}_S{S}_k{k}_{kind}.json').read_text())


def ttt_select(pil, target, capped):
    """Sensitivity variants of the same rule (target, and Goto's cap at T_com when P_S > 0.99)."""
    rows = []
    for k in range(len(P.SB_DTS)):
        best = None
        for S in S_LIST:
            p = float((np.array(pil[(S, k)]['cuts']) >= target).mean())
            v = (S if p > 0.99 else P.mcs99(S, p)) if capped else P.mcs99(S, p)
            if best is None or (v, S) < best:
                best = (v, S)
        mc = float(np.mean(pil[(best[1], k)]['cuts']))
        rows.append(dict(k=k, dt=P.SB_DTS[k], min_ttt=best[0], at_S=best[1], mean_cut_at_S=mc))
    return min(rows, key=lambda r: (r['min_ttt'], -r['mean_cut_at_S'], r['k']))['k'], rows


def select():
    out = {}
    for name in METHODS:
        pil = {(S, k): load(name, S, k, 'pilot') for S in S_LIST for k in range(len(P.SB_DTS))}
        k, rows = RP.sb_select(pil)
        k_same, _ = ttt_select(pil, 33000, False)
        assert k_same == k, 'sensitivity implementation disagrees with run_pub.sb_select at the primary target'
        k33004, rows33004 = ttt_select(pil, 33004, False)
        kcap, rowscap = ttt_select(pil, 33000, True)
        out[name] = dict(selected_k=k, selected_dt=P.SB_DTS[k], rows=rows,
                         sensitivity=dict(target_33004=dict(k=k33004, dt=P.SB_DTS[k33004], rows=rows33004),
                                          goto_capped_33000=dict(k=kcap, dt=P.SB_DTS[kcap], rows=rowscap)))
        print(name, 'selected dt', P.SB_DTS[k], '| target 33004 ->', P.SB_DTS[k33004], '| capped ->', P.SB_DTS[kcap])
        for r in rows:
            print(f"   dt {r['dt']:<4g} min MCS99 {r['min_mcs99']:9.1f} at S = {r['at_S']:<4d} pilot mean cut {r['mean_cut_at_S']:.1f}")
    (OUT / 'selection_A6.json').write_text(json.dumps(out, indent=1, default=float) + '\n')
    return out


def assemble():
    sel = json.loads((OUT / 'selection_A6.json').read_text())
    comp = {}; lines = ['A6 (bSB, dSB: one Delta t for K2000 by time-to-target; every other entry = A2): mean cut / p33000 / MCS99']
    data = {}
    for S in S_LIST:
        a2 = json.loads((A2 / f'S{S}.json').read_text())
        new = {}
        for name, entry in a2.items():
            if name in METHODS:
                k = sel[name]['selected_k']
                new[name] = dict(final=load(name, S, k, 'final'),
                                 setting=f'Goto 2021, dt {P.SB_DTS[k]} (one dt for K2000 by time-to-target, A6)')
            else:
                new[name] = entry
        for name in a2:
            if name not in METHODS:
                assert json.dumps(new[name]) == json.dumps(a2[name]), name
        assert list(new) == list(a2)
        OUT.mkdir(exist_ok=True)
        (OUT / f'S{S}.json').write_text(json.dumps(new) + '\n')
        data[S] = (a2, new)
        lines.append(f'--- S = {S}')
        for n in sorted(new, key=lambda n: -new[n]['final']['mean_cut']):
            f = new[n]['final']
            lines.append(f"{n:20s} mean {f['mean_cut']:9.1f}  p {f['p33000']:.3f}  MCS99 {f['mcs99']:9.0f}   [{new[n].get('setting', 'frozen')}]")
    names = list(data[S_LIST[0]][1])
    mins = {n: min((data[S][1][n]['final']['mcs99'], S) for S in S_LIST) for n in names}
    lines.append('Steps to solution (min over S):')
    for n in sorted(names, key=lambda n: mins[n][0]):
        lines.append(f'  {n:20s} {mins[n][0]:9.0f} at S = {mins[n][1]}')
    lines.append('S = 500 success: ' + ', '.join(f"{n} {data[500][1][n]['final']['p33000']:.3f}" for n in names))
    for name in METHODS:
        c = dict(selected_dt=sel[name]['selected_dt'], per_S={})
        for S in S_LIST:
            a2f = data[S][0][name]['final']; nf = data[S][1][name]['final']
            c['per_S'][S] = dict(A2_setting=data[S][0][name]['setting'], A6_setting=data[S][1][name]['setting'],
                                 mean_cut_A2=a2f['mean_cut'], mean_cut_A6=nf['mean_cut'], p33000_A2=a2f['p33000'],
                                 p33000_A6=nf['p33000'], mcs99_A2=a2f['mcs99'], mcs99_A6=nf['mcs99'])
        a2m = min((data[S][0][name]['final']['mcs99'], S) for S in S_LIST)
        c['mcs99_min_A2'] = dict(value=a2m[0], S=a2m[1]); c['mcs99_min_A6'] = dict(value=mins[name][0], S=mins[name][1])
        comp[name] = c
        lines.append(f"{name}: A6 dt {c['selected_dt']} | MCS99 min A2 {a2m[0]:.0f} (S {a2m[1]}) -> A6 {mins[name][0]:.0f} (S {mins[name][1]})")
    (OUT / 'compare_A6_vs_A2.json').write_text(json.dumps(comp, indent=1, default=float) + '\n')
    t = '\n'.join(lines)
    (OUT / 'claims_A6.txt').write_text(t + '\n')
    print(t)


def check():
    """Module hashes (as synced) and the frozen-module self-test of pub_methods."""
    for f in ('common.py', 'pub_methods.py', 'run_pub.py', 'sb_ttt_A6.py'):
        print(f, hashlib.sha256((HERE / f).read_bytes()).hexdigest())
    P.selftest()
    js = jobs(); assert len({tuple(j['seed']) for j in js}) == len(js) == 100
    print('100 jobs, seeds unique; check OK')


if __name__ == '__main__':
    cmd = sys.argv[1] if len(sys.argv) > 1 else 'check'
    if cmd == 'run':
        run(int(sys.argv[2]) if len(sys.argv) > 2 else 50)
    elif cmd == 'select':
        select()
    elif cmd == 'assemble':
        assemble()
    else:
        check()
