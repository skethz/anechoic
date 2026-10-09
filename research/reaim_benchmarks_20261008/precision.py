"""Coupling-precision study (PROTOCOL.md): the four engine rules on MCP G1-G20, GPP (9) and TSP (6) at K in {2,3,4,6,8}
and full precision, with the hardware arithmetic (hwmodel.py, bit-exact to sca_ref_bias.hpp for K <= 8).
  python3 precision.py select      (laptop) -> precision/schedules.json from the held-out finals
  python3 precision.py run [W]     (gpu-host)   -> precision/runs/<inst>_<rule>_K<K>.json
Schedule per (instance, rule): selected configuration at S_exp = argmin_S MCS99 (ties -> smaller S); if MCS99 is
infinite at every S, the S with the highest final mean quality (ties -> larger S). MCP schedules: research/gset_20261007
pre-registered finals. Quantized schedule: T0, T_fin, q, J_v times alpha_K; kappa_eff and lambda_eff recomputed from J^K."""
import json
import math
import multiprocessing as mp
import os
import sys
import time
from pathlib import Path

for _v in ('NUMBA_NUM_THREADS', 'OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS'):
    os.environ[_v] = '1'
sys.dont_write_bytecode = True
HERE = Path(__file__).resolve().parent
os.environ.setdefault('NUMBA_CACHE_DIR', str(HERE / '.numba_cache'))
sys.path.insert(0, str(HERE))
import numpy as np  # noqa: E402

import hwmodel as HW  # noqa: E402
import problems as PB  # noqa: E402

OUT = HERE / 'precision'
RULES = ('SCA', 'TEC', 'Onsager-kT', 'Onsager-online')
KS = (None, 2, 3, 4, 6, 8)
TRIALS = 1024
GPP_INST = (1, 2, 3, 4, 5, 14, 15, 16, 17)
TSP_INST = ('gr17', 'gr21', 'gr24', 'fri26', 'bayg29', 'bays29')
MCP_INST = tuple(range(1, 21))
GSET_RES = HERE.parent / 'gset_20261007' / 'results' / 'budget'
S_GSET = (250, 500, 1000, 2000, 4000)
S_LIST = {'gpp': (256, 512, 1024, 2048, 4096), 'tsp': (512, 1024, 2048, 4096, 8192)}


def pick(recs):
    """recs: list of (S, final dict, cfg). Returns (S, cfg, rule_text)."""
    fin = [(S, f, c) for S, f, c in recs if f is not None]
    finite = [(f['mcs99'], S, c) for S, f, c in fin if math.isfinite(f['mcs99'])]
    if finite:
        v, S, c = min(finite, key=lambda x: (x[0], x[1]))
        return S, c, f'argmin MCS99 ({v:.0f})'
    key = 'mean_quality' if 'mean_quality' in fin[0][1] else 'mean_cut'
    S, f, c = max(fin, key=lambda x: (x[1][key], x[0]))
    return S, c, f'MCS99 infinite at every S: max final {key} ({f[key]:.4f})'


def select():
    sched = []
    for g in MCP_INST:
        for r in RULES:
            recs = []
            for S in S_GSET:
                d = json.loads((GSET_RES / f'G{g}' / f"{r.replace(' ', '_')}_S{S}.json").read_text())
                recs.append((S, d['final'], d['selected_cfg']))
            S, c, why = pick(recs)
            sched.append(dict(problem='mcp', instance=g, rule=r, S=S, cfg=c, why=why, source='research/gset_20261007 pre-registered'))
    for prob, insts in (('gpp', GPP_INST), ('tsp', TSP_INST)):
        for i in insts:
            tag = f'G{i}' if prob == 'gpp' else i
            for r in RULES:
                recs = []
                for S in S_LIST[prob]:
                    d = json.loads((HERE / 'results' / prob / tag / f"{r.replace(' ', '_')}_S{S}.json").read_text())
                    recs.append((S, d['final'], d['selected_cfg']))
                S, c, why = pick(recs)
                sched.append(dict(problem=prob, instance=i, rule=r, S=S, cfg=c, why=why, source='this study held-out finals'))
    OUT.mkdir(exist_ok=True)
    (OUT / 'schedules.json').write_text(json.dumps(sched, indent=1) + '\n')
    print(f'{len(sched)} schedules written')


def problem(prob, inst):
    if prob == 'mcp':
        bkv = json.loads((HERE.parent / 'gset_20261007' / 'bkv.json').read_text())['instances']
        return PB.mcp(inst, bkv[f'G{inst}'])
    if prob == 'gpp':
        ref = json.loads((HERE / 'data' / 'gpp_reference.json').read_text())['instances'][f'G{inst}']
        return PB.gpp(inst, 4, ref=ref['R'], target=ref['target'])
    return PB.tsp(inst)


def rel_params(prob, cfg, P):
    """Relative (sigma_T / factor) parameters of a selected configuration."""
    a = P.sigma_T; st = P.stats(); fam = cfg['family']
    r = dict(family=fam, T0=cfg['T0'] / a, tfin=cfg['tfin'] / a, q=cfg['q'] / a, ramp=bool(cfg.get('ramp', False)))
    if fam == 'tec':
        r['jv'] = cfg['jv'] / a
    if fam == 'tecT':
        r['kappa'] = cfg['kappa'] / st['kappa_factor']
    if fam == 'onsager':
        r['lam'] = cfg['lam'] / st['lam_factor']
    return r


def hw_schedule(rel, Q, S):
    """Schedule in the units of the (quantized) problem Q (hwmodel / mb_common.hpp Sched fields)."""
    a = Q.sigma_T; st = Q.stats()
    sch = dict(t0=rel['T0'] * a, t1=rel['tfin'] * a, S=int(S), q=rel['q'] * a, ramp=rel['ramp'], lam=0.0, jv=0.0, kappa=0.0)
    if rel['family'] == 'tec':
        sch['jv'] = rel['jv'] * a
    if rel['family'] == 'tecT':
        sch['kappa'] = rel['kappa'] * st['kappa_factor']
    if rel['family'] == 'onsager':
        sch['lam'] = rel['lam'] * st['lam_factor']
    return sch


PCODE = {'mcp': 0, 'gpp': 1, 'tsp': 2}


def icode(prob, inst):
    return int(inst) if prob != 'tsp' else 100 + TSP_INST.index(inst)


def job(args):
    k, s, K = args
    tag = f"{s['problem']}_{s['instance']}_{s['rule'].replace(' ', '_')}_K{K or 'full'}"
    path = OUT / 'runs' / f'{tag}.json'
    if path.exists():
        return None
    t0 = time.time()
    P = problem(s['problem'], s['instance'])
    rel = rel_params(s['problem'], s['cfg'], P)
    Q, alpha = P.quantize(K)
    sch = hw_schedule(rel, Q, s['S'])
    seed = 20261008 + 1000 * PCODE[s['problem']] + 10 * (icode(s['problem'], s['instance']) % 100) + RULES.index(s['rule'])
    out = HW.run(Q, sch, seed, np.arange(TRIALS))
    ev = P.evaluate(out['spins'].astype(np.float32))
    rec = dict(problem=s['problem'], instance=P.name, rule=s['rule'], K=K, alpha=alpha, S=s['S'], schedule=sch, rel=rel,
               seed=seed, trials=TRIALS, checks=HW.table_checks(sch, K), qstats=Q.stats(),
               quality=ev['quality'].tolist(), feasible=ev['feasible'].tolist(), success=ev['success'].tolist(),
               value=ev['value'].tolist(), flips=out['flips'].tolist(), sec=time.time() - t0,
               mean_quality=float(ev['quality'].mean()), p_feasible=float(ev['feasible'].mean()),
               p_target=float(ev['success'].mean()))
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_suffix('.tmp'); tmp.write_text(json.dumps(rec) + '\n'); tmp.rename(path)
    return (f"{P.name:7s} {s['rule']:15s} K={str(K):4s} S={s['S']:5d} q={rec['mean_quality']:.4f} feas={rec['p_feasible']:.3f} "
            f"p={rec['p_target']:.3f} ({rec['sec']:.0f}s)")


def run(W):
    sched = json.loads((OUT / 'schedules.json').read_text())
    jobs = [(k, s, K) for k, s in enumerate(sched) for K in KS]
    jobs.sort(key=lambda j: -(j[1]['S'] * (900 if j[1]['problem'] == 'tsp' else 800)))
    print(f'precision: {len(jobs)} jobs, {W} workers', flush=True)
    t0 = time.time()
    with mp.Pool(W) as pool:
        for n, msg in enumerate(pool.imap_unordered(job, jobs, chunksize=1)):
            if msg:
                print(f'[{n + 1}/{len(jobs)} {time.time() - t0:.0f}s] {msg}', flush=True)
    print(f'precision done {time.time() - t0:.0f}s', flush=True)


if __name__ == '__main__':
    if sys.argv[1] == 'select':
        select()
    else:
        run(int(sys.argv[2]) if len(sys.argv) > 2 else 48)
