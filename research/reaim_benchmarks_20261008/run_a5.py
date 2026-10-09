"""Addendum 5 (the authors' decisions after the fairness audit, research/fairness_audit_20261008): Table 8b re-runs.
  ReAIM ASA with ReAIM Table I's per-problem k ("Best flips k per iteration": MCP >= 16, GPP 6, TSP 1):
     MCP k in {128, 256, 512, 1024} capped at N (the >= 16 set with which the audit reproduces ReAIM's K2000 Table VII),
     GPP {6}, TSP {1}; Table II T, F per Sec. IV-C, ITER_trial/run 32/96, x_best output (Algorithm 3), noise-free.
     finals: 256 runs, SeedSequence([20261008, 80, p, i, m, s]) -> results_a5/final/reaim_t1/<prob>/<inst>/ReAIM_ASA_S<S>.json
  bSB, dSB: one dt per instance by time-to-target (Goto 2021: dt "set to the best value for each problem among five
     values (0.25, 0.5, 0.75, 1, and 1.25)"; "N_step is also optimized for TTT"; TTT target 99 % of the best known
     = our MCS99 target). For every dt in that set and every S of the grid:
       pilot 256 runs, SeedSequence([20261008, 82, p, i, m, s, k]) -> results_a5/pilot/sb_ttt/<prob>/<inst>/<m>_dt<dt>_S<S>.json
       final 256 runs, SeedSequence([20261008, 83, p, i, m, s, k]) -> results_a5/final/sb_ttt/<prob>/<inst>/<m>_dt<dt>_S<S>.json
     selection ('select'): dt minimising min over S of the pilot MCS99; ties (also all infinite) to the higher pilot mean
     quality at the S attaining the minimum (the largest S when all are infinite), then to the smaller dt
     -> results_a5/sb_ttt_selection.json. TSP: ancillary spin (Addendum 4.1, ancilla.py).
Usage: python3 run_a5.py run [W] | select | check"""
import os

for _v in ('NUMBA_NUM_THREADS', 'OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'VECLIB_MAXIMUM_THREADS', 'MKL_NUM_THREADS'):
    os.environ[_v] = '1'

import json  # noqa: E402
import math  # noqa: E402
import sys  # noqa: E402
import time  # noqa: E402
from pathlib import Path  # noqa: E402

sys.dont_write_bytecode = True
HERE = Path(__file__).resolve().parent
os.environ.setdefault('NUMBA_CACHE_DIR', str(HERE / '.numba_cache'))
sys.path.insert(0, str(HERE))
import numpy as np  # noqa: E402

import ancilla as AN  # noqa: E402
import grids_a2 as GA  # noqa: E402
import run_a2 as RA  # noqa: E402
import run_a3 as R3  # noqa: E402
import settings_a3 as ST  # noqa: E402
import solvers as SV  # noqa: E402
import solvers_a3b as SV3B  # noqa: E402

OUT = HERE / 'results_a5'
SEED = 20261008
B = 256
REAIM_KSET_T1 = {'mcp': (128, 256, 512, 1024), 'gpp': (6,), 'tsp': (1,)}


def reaim_cfg(prob, P, S):
    T0, T1 = ST.REAIM_TABLE_II[prob]
    ks = tuple(sorted(set(int(min(P.N, k)) for k in REAIM_KSET_T1[prob])))
    return dict(family='ReAIM', kset=ks, T0=T0, T1=T1, F=ST.REAIM_F[prob], S=S,
                note='ReAIM Table I per-problem k (MCP >=16: {128,256,512,1024} capped at N; GPP 6; TSP 1), Table II T, '
                     'F per Sec. IV-C, ITER_trial/run 32/96 (not given), x_best output, noise-free')


def save(path, rec):
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_suffix('.tmp'); tmp.write_text(json.dumps(rec, default=list) + '\n'); tmp.rename(path)


def reaim_job(prob, inst, S):
    path = OUT / 'final' / 'reaim_t1' / prob / RA.tag(prob, inst) / f'ReAIM_ASA_S{S}.json'
    if path.exists():
        return None
    t0 = time.time()
    P = RA.problem(prob, inst)
    mi = R3.METHODS.index('ReAIM ASA'); si = GA.S_LIST[prob].index(S); pc = RA.PCODE[prob]; ic = RA.icode(prob, inst)
    cfg = reaim_cfg(prob, P, S)
    s, _, _ = SV3B.run_method(P, cfg, B, np.random.SeedSequence([SEED, 80, pc, ic, mi, si]))
    ev = P.evaluate(s)
    fin, q, feas = R3.stats(P, prob, ev, None, None, S)
    rec = dict(problem=prob, instance=P.name, method='ReAIM ASA', S=S, cfg=cfg, final=fin, quality=q.tolist(),
               values=ev['value'].tolist(), feasible_runs=feas.tolist(), sec=time.time() - t0)
    if prob == 'gpp':
        rec['imbalance'] = ev['imbalance'].tolist()
    save(path, rec)
    return f"ReAIM-T1 {prob} {RA.tag(prob, inst):6s} S={S:<5d} q={fin['mean_quality']:.4f} feas={fin['p_feasible']:.3f} p={fin['p_target']:.3f} ({rec['sec']:.0f}s)"


def sb_job(kind, method, prob, inst, S, k):
    dt = ST.SB_DT_SET[k]
    path = OUT / kind / 'sb_ttt' / prob / RA.tag(prob, inst) / f'{method}_dt{dt:g}_S{S}.json'
    if path.exists():
        return None
    t0 = time.time()
    P = RA.problem(prob, inst)
    Q = AN.with_ancilla(P) if prob == 'tsp' else P
    mi = R3.METHODS.index(method); si = GA.S_LIST[prob].index(S); pc = RA.PCODE[prob]; ic = RA.icode(prob, inst)
    cfg = dict(family=method, dt=dt, xi=1.0, S=S,
               note='Goto 2021: a0=1, c0=0.5/(<J> sqrt N), dt from {0.25,...,1.25} by TTT per instance'
                    + ('; ancillary spin for the fields (A4.1)' if prob == 'tsp' else ''))
    rng = np.random.default_rng(np.random.SeedSequence([SEED, 82 if kind == 'pilot' else 83, pc, ic, mi, si, k]))
    s = SV.sb(Q, B, cfg, rng)
    if prob == 'tsp':
        s = AN.decode(s)
    ev = P.evaluate(s)
    fin, q, feas = R3.stats(P, prob, ev, None, None, S)
    rec = dict(kind=kind, problem=prob, instance=P.name, method=method, S=S, dt=dt, cfg=cfg, final=fin,
               quality=q.tolist(), values=ev['value'].tolist(), feasible_runs=feas.tolist(), sec=time.time() - t0)
    save(path, rec)
    return (f"SB-TTT {kind:5s} {method} {prob} {RA.tag(prob, inst):6s} dt={dt:<4g} S={S:<5d} q={fin['mean_quality']:.4f} "
            f"p={fin['p_target']:.3f} ({rec['sec']:.0f}s)")


def job(args):
    if args[0] == 'reaim':
        return reaim_job(*args[1:])
    return sb_job(*args)


def jobs():
    J = []
    for prob, insts in GA.INSTANCES.items():
        for inst in insts:
            for S in GA.S_LIST[prob]:
                J.append(('reaim', prob, inst, S))
                for method in ('bSB', 'dSB'):
                    for k in range(len(ST.SB_DT_SET)):
                        J.append(('pilot', method, prob, inst, S, k))
                        J.append(('final', method, prob, inst, S, k))
    size = {'mcp': 800, 'gpp': 800, 'tsp': 900}
    J.sort(key=lambda j: -(size[j[1] if j[0] == 'reaim' else j[2]] * (j[3] if j[0] == 'reaim' else j[4])
                           * (4 if j[0] == 'reaim' else 1) * (3 if (j[0] != 'reaim' and j[2] == 'tsp') else 1)))
    return J


def select():
    out = {}
    for prob, insts in GA.INSTANCES.items():
        for inst in insts:
            for method in ('bSB', 'dSB'):
                rows = []
                for k, dt in enumerate(ST.SB_DT_SET):
                    m = []
                    for S in GA.S_LIST[prob]:
                        d = json.loads((OUT / 'pilot' / 'sb_ttt' / prob / RA.tag(prob, inst) / f'{method}_dt{dt:g}_S{S}.json').read_text())
                        m.append((float(d['final']['mcs99']), S, d['final']['mean_quality']))
                    fin = [x for x in m if math.isfinite(x[0])]
                    if fin:
                        best = min(fin, key=lambda x: (x[0], x[1]))
                    else:
                        best = (math.inf, m[-1][1], m[-1][2])
                    rows.append(dict(dt=dt, ttt=best[0], S_at_min=best[1], q_at_S=best[2],
                                     pilot_mcs99={str(x[1]): x[0] for x in m}, pilot_quality={str(x[1]): x[2] for x in m}))
                sel = min(range(len(rows)), key=lambda k: (rows[k]['ttt'], -rows[k]['q_at_S'], k))
                out[f'{prob}/{RA.tag(prob, inst)}/{method}'] = dict(selected_dt=rows[sel]['dt'], rows=rows)
    (OUT / 'sb_ttt_selection.json').write_text(json.dumps(out, indent=1, default=float) + '\n')
    for key, v in out.items():
        print(key, 'dt =', v['selected_dt'], 'TTT =', [f"{r['ttt']:.0f}" if math.isfinite(r['ttt']) else 'inf' for r in v['rows']])
    return out


def check():
    """Configuration plumbing: the per-problem k sets as passed to reaim_best, and a short run of each."""
    ok = True
    for prob, inst in (('mcp', 1), ('gpp', 1), ('tsp', 'gr17')):
        P = RA.problem(prob, inst)
        cfg = reaim_cfg(prob, P, 64)
        exp = {'mcp': (128, 256, 512, 800), 'gpp': (6,), 'tsp': (1,)}[prob]
        good = cfg['kset'] == exp
        s, _, _ = SV3B.run_method(P, cfg, 4, np.random.SeedSequence([1]))
        good &= s.shape == (4, P.N)
        print(prob, inst, 'kset', cfg['kset'], 'expected', exp, 'OK' if good else 'FAIL', flush=True)
        ok &= bool(good)
    print('ALL OK' if ok else 'FAIL')
    return ok


if __name__ == '__main__':
    cmd = sys.argv[1]
    if cmd == 'check':
        sys.exit(0 if check() else 1)
    if cmd == 'select':
        select(); sys.exit(0)
    W = int(sys.argv[2]) if len(sys.argv) > 2 else 64
    RA.run_pool(jobs(), job, W, 'A5 Table 8b re-runs (ReAIM Table I k; bSB/dSB TTT dt)')
