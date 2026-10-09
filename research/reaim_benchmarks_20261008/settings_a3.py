"""Addendum 3: the eight baselines at the settings their own papers describe (no tuning by us), and the two Onsager
forms at their unchanged A1 selections. Writes settings_a3.json (every instance's resolved setting, with provenance).

Sources (local copies in sources/papers/, SHA-256 in sources/papers/pdf_sha256.txt):
  STATICA  Yamamoto et al., IEEE JSSC 56(1):165-178, 2021 (yamamoto2021statica): Alg. 1, Eq. (7), Sec. V-C, Table I,
           Fig. 24, Fig. 25.
  TEC      Du et al., arXiv:2608.21753, 2026 (du2026tec): Eqs. (2)-(4), Fig. 3 and its text.
  APC-SCA  Okonogi et al., IEICE Trans. Inf. & Syst. E106-D(12):1969-1978, 2023 (okonogi2023apc): Alg. 2, Eq. (3),
           Table 2, Sec. 4.2, Sec. 4.4 / Table 5.
  ReAIM    Chiang et al., ISCA 2024 (10609617): Algs. 2-3, Steps 0-5 (Sec. IV-B), Sec. IV-C, Table I, Table II, Fig. 3.
  aSB      Goto, Tatsumura, Dixon, Sci. Adv. 5:eaav2372, 2019 (goto2019combinatorial): Eqs. (14)-(17), Fig. 2 caption,
           Methods ("we set dt = 0.9 and M = 2 in Fig. 2").
  bSB/dSB  Goto et al., Sci. Adv. 7:eabe7953, 2021 (goto2021high): Methods "Parameter setting".
  SA       D-Wave Neal (ReAIM's CPU baseline [15]): dwave-samplers 1.2.0 SimulatedAnnealingSampler (the implementation
           behind `neal` 0.6.0; sources/papers/dwave_samplers_1.2.0_sa_sampler.py, cpu_sa.cpp): default beta range
           `_default_ising_beta_range`, geometric beta schedule, sequential sweeps in variable order, Metropolis.
Transfer rule for settings a paper gives only for K2000 (STATICA, TEC): field-valued values are multiplied by
sigma_inst / sigma_K2000, sigma = sqrt(mean_i sum_j J_ij^2) (the RMS local field; K2000: sqrt(1999)), the rule of the
project's pre-registered G-set study. This is ours and is flagged as such.
"""
import json
import math
import sys
from pathlib import Path

import numpy as np

sys.dont_write_bytecode = True
HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))

SIGMA_K2000 = math.sqrt(1999.0)
REAIM_TABLE_II = {'mcp': (1.0, 0.1), 'gpp': (1.0, 0.01), 'tsp': (0.5, 0.1)}   # T_init / T_final
REAIM_F = {'mcp': 'max', 'gpp': 'max', 'tsp': 'min'}                           # Sec. IV-C (outcome of Step 4)
REAIM_KSET = (1, 2, 6, 16)       # every k ReAIM reports as best/used: Table I (>=16, 6, 1), Fig. 3 (16, 6, 2)
SB_DT_SET = (0.25, 0.5, 0.75, 1.0, 1.25)                                       # Goto 2021, Methods
BASELINES = ('SA', 'SCA', 'TEC', 'APC-SCA', 'ReAIM ASA', 'aSB', 'bSB', 'dSB')
OURS = ('Onsager-kT', 'Onsager-online')


def instance_constants(P):
    """lambda_max(-J), SD of the elements of J, min/max nonzero |J|, RMS field, Neal's default beta range."""
    J = P.dense(); N = P.N; b = P.b.astype(np.float64)
    lam = float(np.linalg.eigvalsh(-J)[-1])
    sd = float(np.sqrt((J * J).mean() - J.mean() ** 2))                   # all N^2 elements (diagonal zeros included)
    nz = np.abs(J[J != 0])
    minJ, maxJ = float(nz.min()), float(nz.max())
    sigma_full = math.sqrt(float((J * J).sum(1).mean()))
    absJ = np.abs(J)
    sum_abs = np.abs(b) + absJ.sum(1)
    big = np.where(absJ > 0, absJ, np.inf).min(1)
    min_abs = np.where(np.abs(b) > 0, np.minimum(np.abs(b), big), big)
    max_eff = float(sum_abs.max()); min_eff = float(min_abs.min()); n_min = int(np.sum(min_abs == min_eff))
    hot = math.log(2) / (2 * max_eff); cold = math.log(n_min / 0.01) / (2 * min_eff)
    return dict(N=N, lam_max_minusJ=lam, sd_J=sd, rms_J=P.stats()['sigmaJ_goto'], minJ=minJ, maxJ=maxJ,
                sigma_full=sigma_full, neal=dict(max_effective_field=max_eff, min_effective_field=min_eff,
                                                 number_min_gaps=n_min, hot_beta=hot, cold_beta=cold))


def baseline_cfgs(prob, P, method, S, C):
    """List of configurations (one, or the five dt candidates of bSB/dSB) with provenance notes."""
    r = C['sigma_full'] / SIGMA_K2000
    N = P.N
    if method == 'SA':
        return [dict(family='SA', T0=1.0 / C['neal']['hot_beta'], T1=1.0 / C['neal']['cold_beta'], S=S,
                     note='Neal default beta range (dwave-samplers 1.2.0), geometric, S sweeps')]
    if method == 'SCA':
        return [dict(family='plain', q=4.0 * r, T0=40.0 * r, tfin=5.0 * r, S=S,
                     note='STATICA K2000 setting q=4 (constant), T 40->5 (Fig. 25, Table II), x sigma/sigma_K2000')]
    if method == 'TEC':
        return [dict(family='tec_sig', jv=30.0 * r, q=0.0, T0=100.0 * r, tfin=0.1 * r, cb=2.0, S=S,
                     note='TEC K2000: Jv = 30J optimum (Fig. 3b), T 100J->0.1J (Fig. 3 text), Glauber Eq. (4), '
                          'x sigma/sigma_K2000; geometric schedule (shape not stated)')]
    if method == 'APC-SCA':
        lam2 = C['lam_max_minusJ'] / 2
        if prob == 'mcp':
            return [dict(family='apc_sig', q_reset=lam2, r_q=0.45, q_lim=0.0, T0=10.0, tfin=0.1, cb=1.0, S=S,
                         note='APC max-cut: T 10->0.1 (Table 2), r_q=0.45, q_limit=0 (Sec. 4.2), q_i(1)=reset=lambda/2 (Alg. 2)')]
        return [dict(family='apc_sig', q_reset=lam2, r_q=0.8, q_lim=0.0, T0=0.2 * N * C['minJ'], tfin=0.1 * C['maxJ'],
                     cb=1.0, S=S, note=('APC TSP: T 0.2 N min|J| -> 0.1 max|J| (Table 2), r_q=0.8 (Sec. 4.4, Table 5), '
                                        'q_limit=0, q=lambda/2' + (' [GPP not covered: TSP column used]' if prob == 'gpp' else '')))]
    if method == 'ReAIM ASA':
        T0, T1 = REAIM_TABLE_II[prob]
        return [dict(family='ReAIM', kset=tuple(int(min(N, k)) for k in REAIM_KSET), T0=T0, T1=T1, F=REAIM_F[prob], S=S,
                     note='ReAIM Table II T_init/T_final, F per Sec. IV-C, k candidates {1,2,6,16} (paper values), '
                          'ITER_trial/run 32/96 (not given), noise-free')]
    if method == 'aSB':
        return [dict(family='aSB', dt=0.9, M=2, xi=C['rms_J'] / C['sd_J'], S=S,
                     note='Goto 2019: K=D=1, xi0=0.7/(SD(J) sqrt N), p linear 0->1, dt=0.9, M=2 (Fig. 2, Methods)')]
    if method in ('bSB', 'dSB'):
        return [dict(family=method, dt=dt, xi=1.0, S=S,
                     note='Goto 2021: a0=1, c0=0.5/(<J> sqrt N), a linear 0->a0, dt best of {0.25,...,1.25} (Methods)')
                for dt in SB_DT_SET]
    raise ValueError(method)


def main():
    import run_a2 as RA                       # problem constructors only
    import grids_a2 as GA
    out = {}
    for prob, insts in GA.INSTANCES.items():
        for inst in insts:
            P = RA.problem(prob, inst)
            C = instance_constants(P)
            key = f'{prob}/{RA.tag(prob, inst)}'
            out[key] = dict(constants=C, settings={m: baseline_cfgs(prob, P, m, GA.S_LIST[prob][-1], C) for m in BASELINES})
            s = out[key]['settings']
            print(f"{key:12s} N={P.N:4d} lam/2={C['lam_max_minusJ'] / 2:9.1f} sd={C['sd_J']:.3f} rms={C['rms_J']:.3f} "
                  f"sig={C['sigma_full']:8.1f} Neal T {s['SA'][0]['T0']:.4g}->{s['SA'][0]['T1']:.4g} | SCA q={s['SCA'][0]['q']:.3g} "
                  f"T {s['SCA'][0]['T0']:.4g}->{s['SCA'][0]['tfin']:.4g} | TEC Jv={s['TEC'][0]['jv']:.3g} | APC T "
                  f"{s['APC-SCA'][0]['T0']:.4g}->{s['APC-SCA'][0]['tfin']:.4g}", flush=True)
    (HERE / 'settings_a3.json').write_text(json.dumps(out, indent=1) + '\n')


if __name__ == '__main__':
    main()
