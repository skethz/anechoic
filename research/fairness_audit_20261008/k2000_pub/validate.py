"""PROTOCOL_PUB_V1.md: does each baseline implementation, at its paper's settings, reproduce a number its own paper reports?
Fresh seeds SeedSequence([20261008, 601, case index]). Usage: python3 validate.py run [workers] | report
Cases:
  statica_long   K2000, clipped Eq. (7), q 4, T 40 -> 5, S 1560, final state, 1,024 runs  vs JSSC 2021 Table II / Fig. 22:
                 P_a = 0.77 (100 runs; exact 95% CI 0.675-0.848), mean cut 33,073
  statica_short  same with T 30 -> 5, S 560                                                vs P_a = 0.07 (0.029-0.139), mean cut 32,750
  reaim_4096     K2000, ASA k {128,256,512,1024}, F = max, T 1 -> 0.1, ITER 32/96, output x_best (Alg. 3), 512 runs
                 vs ISCA 2024 Table VII: P_a = 0.47 at t_a = 0.15 ms
  reaim_6400     same, 6,400 iterations (6,400/4,096 = 1.56 ~ 0.23/0.15 = 1.53)          vs P_a = 0.8 at t_a = 0.23 ms
  reaim_cop_6400 the COP study's k set {1, 2, 6, 16}, otherwise as reaim_6400, 256 runs   vs P_a = 0.8
  asb_186        K2000, aSB Delta t 0.9, M 2, xi0 0.7/sqrt(N), N_step 186, per-run divergence = failure, 1,024 runs
                 vs Goto et al. 2019 FPGA K2000 point quoted in JSSC 2021 Table II / ISCA 2024 Table VII: 0.5 ms, P_a = 0.04
  apc_table3_G   G22, G30, G32, G35: the APC paper's 'fine-tuned SCA' (Algorithm 1, logistic, q(s) exponential from q_init
                 to q_final of its Table 3, T 10 -> 0.1 exponential, S 1000, output argmin_s H), 128 rounds
                 vs APC Table 3 obtained energies -6545.5, -6638.9, -2726.3, -3430.1
  apc_G22        G22, APC-SCA Algorithm 2 with r_q 0.45, q_limit 0, T 10 -> 0.1, S 1000, output argmin, 128 rounds
                 vs the paper's claim (Sect. 4.2, Figs. 3-4): lower average energy than fine-tuned SCA (-6545.5)
"""
import json
import math
import os
import sys
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
import pub_methods as P  # noqa: E402
from common import M, ROOT  # noqa: E402

HERE = Path(__file__).resolve().parent
OUT = HERE / 'results_validation'
SEED = 20261008
APC_TABLE3 = {'G22': (31, 1, -6545.5), 'G30': (2, 2, -6638.9), 'G32': (26, 3, -2726.3), 'G35': (27, 1, -3430.1)}
CASES = ['statica_long', 'statica_short', 'reaim_4096', 'reaim_6400', 'reaim_cop_6400', 'asb_186',
         'apc_table3_G22', 'apc_table3_G30', 'apc_table3_G32', 'apc_table3_G35', 'apc_G22']


def gset_J(name):
    p = ROOT / 'research/gset_20261007/data' / name
    lines = p.read_text().split('\n'); n, m = map(int, lines[0].split())
    a = np.array([ln.split() for ln in lines[1:] if ln.strip()], dtype=np.int64)
    J = np.zeros((n, n), np.float32)
    J[a[:, 0] - 1, a[:, 1] - 1] = -a[:, 2]; J[a[:, 1] - 1, a[:, 0] - 1] = -a[:, 2]
    return J


def logistic_sca_best(J, B, S, q_sched, T0, T1, rng, apc=None):
    """APC paper Algorithm 1 (q_sched per step) or Algorithm 2 (apc = (lam, r_q, q_limit)); output argmin_s H incl. s = 1."""
    N = len(J); T = P.geom_sched(T0, T1, S)
    s = rng.choice(np.array([-1.0, 1.0], dtype=np.float32), size=(B, N)); h = s @ J
    H = P.energy_from_field(s, h); best = H.copy()
    q = np.full((B, N), apc[0] / 2 if apc else 0.0)
    for t in range(S):
        qq = q if apc else q_sched[t]
        z = (s * h).astype(np.float64) + qq
        pf = 1.0 / (1.0 + np.exp(np.clip(z / T[t], -700, 700)))
        flip = pf > rng.random(s.shape, dtype=np.float32)
        if apc:
            q = np.where(flip, apc[0] / 2, np.maximum(q * apc[1], apc[2]))
        d = np.where(flip, -2 * s, 0).astype(np.float32)
        s += d; h += d @ J
        best = np.minimum(best, P.energy_from_field(s, h))
    return best


def reaim_best(J, B, iters, kset, rng):
    """pilot.asa (Algorithms 2-3, F = max, T 1 -> 0.1, ITER 32/96); returns (H_best after run phases, final H)."""
    best, final = M.reaim.asa(J, B, iters, kset, 32, 96, rng, 1.0, 0.1, np.max)
    return best, final


def job(ci):
    np.seterr(all='ignore')
    case = CASES[ci]
    rng = np.random.default_rng(np.random.SeedSequence([SEED, 601, ci]))
    out = dict(case=case)
    if case.startswith('statica'):
        J, sumw = M.m.load()
        S, T0 = (1560, 40.0) if case == 'statica_long' else (560, 30.0)
        B = 1024
        s0 = rng.choice(np.array([-1.0, 1.0], dtype=np.float32), size=(B, len(J)))
        s, fl = P.abl.run(J, s0, dict(family='plain', q=4.0, T0=T0, S=S), rng)
        c = P.m.cut(J, sumw, s)
        out.update(runs=B, k=int((c >= 33000).sum()), mean_cut=float(c.mean()), sd_cut=float(c.std(ddof=1)))
    elif case.startswith('reaim'):
        J, sumw = M.m.load()
        iters = 4096 if case == 'reaim_4096' else 6400
        kset = (1, 2, 6, 16) if case == 'reaim_cop_6400' else (128, 256, 512, 1024)
        B = 256 if case == 'reaim_cop_6400' else 512
        best, final = reaim_best(J, B, iters, kset, rng)
        cb = (sumw - best) / 2; cf = (sumw - final) / 2
        out.update(runs=B, kset=kset, iters=iters, k=int((cb >= 33000).sum()), k_final=int((cf >= 33000).sum()),
                   mean_cut_best=float(cb.mean()))
    elif case == 'asb_186':
        J, sumw = M.m.load()
        import asb_perrun as A
        B = 1024
        s, E, ok = A.asb_traj_flags(J, B, dict(family='aSB', dt=0.9, xi=1.0, M=2, S=186), rng)
        c = (sumw - E[-1]) / 2
        out.update(runs=B, diverged=int((~ok).sum()), k=int(((c >= 33000) & ok).sum()),
                   mean_cut_nondiverged=float(c[ok].mean()))
    elif case.startswith('apc_table3'):
        g = case.split('_')[-1]; J = gset_J(g)
        lam = float(np.linalg.eigvalsh(-J.astype(np.float64))[-1])
        a, b, ref = APC_TABLE3[g]
        q0, q1 = a / 64 * lam, b / 64 * lam
        S = 1000; qs = q0 * (q1 / q0) ** (np.arange(S) / (S - 1))
        best = logistic_sca_best(J, 128, S, qs, 10.0, 0.1, rng)
        out.update(runs=128, lam=lam, q_init=q0, q_final=q1, mean_best_H=float(best.mean()), sd_best_H=float(best.std(ddof=1)),
                   paper=ref)
    elif case == 'apc_G22':
        J = gset_J('G22'); lam = float(np.linalg.eigvalsh(-J.astype(np.float64))[-1])
        best = logistic_sca_best(J, 128, 1000, None, 10.0, 0.1, rng, apc=(lam, 0.45, 0.0))
        out.update(runs=128, lam=lam, mean_best_H=float(best.mean()), sd_best_H=float(best.std(ddof=1)), fine_tuned_sca_paper=-6545.5)
    return out


def run(workers):
    OUT.mkdir(exist_ok=True)
    with ProcessPoolExecutor(max_workers=workers) as ex:
        for r in ex.map(job, range(len(CASES))):
            (OUT / f"{r['case']}.json").write_text(json.dumps(r, indent=1) + '\n'); print(json.dumps(r), flush=True)
    report()


def report():
    L = []
    for case in CASES:
        f = OUT / f'{case}.json'
        if not f.exists():
            continue
        r = json.loads(f.read_text())
        if 'k' in r:
            p = r['k'] / r['runs']; lo, hi = P.wilson(r['k'], r['runs'])
            r['p'] = p; r['wilson95'] = [lo, hi]
        L.append(json.dumps(r))
    t = '\n'.join(L)
    print(t)
    (OUT / 'report.txt').write_text(t + '\n')


if __name__ == '__main__':
    if sys.argv[1] == 'run':
        run(int(sys.argv[2]) if len(sys.argv) > 2 else 11)
    else:
        report()
