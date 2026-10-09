"""Baselines at their PUBLISHED settings for the K2000 algorithm comparison (PROTOCOL_PUB.md). No baseline is tuned by us.

Settings (sources in research/fairness_audit_20261008/sources/; see PROTOCOL_PUB.md for the quoted text):
  SA (Neal)       dwave-samplers 1.5.0 / dwave-neal 0.6.0 SimulatedAnnealingSampler defaults: beta range from
                  _default_ising_beta_range(h, J) = [ln2/(2*1999), ln(2000/0.01)/2] on K2000, geometric beta schedule
                  with one beta per sweep (num_sweeps_per_beta = 1), Metropolis, sequential order, random start,
                  num_sweeps = S. Kernel: methods._sa_kernel (sequential Metropolis sweeps in index order).
  SCA (STATICA)   Yamamoto et al., JSSC 2021: clipped Eq. (7), q = 4, T_fin = 5, geometric T, T_init from STATICA's
                  K2000 grid {30, 40, 50}, chosen per budget by STATICA's figure of merit (TTS), on pilot runs.
  TEC (Du et al.) arXiv:2608.21753: simulated annealing of a p-bit network with Glauber (heat-bath) updates and the
                  temporal field J_v*sigma_i(t-1) (Eq. 4); K2000 settings of Fig. 3: k_B T from 100 J to 0.1 J,
                  J_v = +30 J (the paper's optimum on K2000). Not stated, declared here: one annealing cycle = one
                  sequential sweep (index order), sigma(t-1) = configuration at the end of the previous cycle (the
                  random start for the first cycle), geometric temperature schedule.
  TEC on STATICA  supplementary: STATICA's SCA (as above) plus TEC's published J_v = +30 (field h + J_v s(t-1)).
  APC-SCA         Okonogi et al., IEICE Trans. E106-D(12), 2023, Algorithm 2: logistic flip probability
                  sigmoid(-(h_i s_i + q_i)/T), q_i(1) = lambda/2 (lambda = largest eigenvalue of -J), q_i reset to
                  lambda/2 after a flip, else q_i <- max(q_i r_q, q_limit); max-cut settings: r_q = 0.45, q_limit = 0
                  (Sect. 4.2), exponential T from 10 to 0.1 (Table 2, Eq. 3).
  ReAIM ASA       Chiang et al., ISCA 2024, Algorithms 2-3, noise-free, Max-Cut settings T 1 -> 0.1, F = max.
                  Unreported for K2000 (declared): k set {128, 256, 512, 1024}, ITER_trial 32, ITER_run 96, FIFO
                  initialised with N: the values with which research/reaim_reproduction_20261003 reproduces ReAIM's
                  published K2000 success probabilities (0.47 at about 4,100 iterations, 0.80 at about 6,400).
  aSB             Goto et al., Sci. Adv. 2019: K = Delta = 1, xi0 = 0.7 Delta/(sigma sqrt N) (sigma = 1 on K2000), p(t)
                  linear 0 -> 1, modified symplectic Euler with Delta t = 0.9 and M = 2 (the K2000 setting of Fig. 2),
                  x(0) = 0, y(0) uniform in (-0.1, 0.1). Integrator: methods.asb_traj.
  bSB, dSB        Goto et al., Sci. Adv. 2021: a0 = 1, a(t) linear 0 -> a0, c0 = 0.5/(<J> sqrt N), inelastic walls,
                  Delta t = the best of {0.25, 0.5, 0.75, 1, 1.25} for the problem, with N_step optimised for the
                  target (here: minimum over S of the steps to solution on pilot runs). Integrator: methods.sb_traj.
"""
import math
import sys
from pathlib import Path

import numpy as np
import numba as nb

sys.path.insert(0, str(Path(__file__).resolve().parent))
from common import M, ROOT  # noqa: E402

abl = M.abl
m = M.m

NEAL_BETA = (0.00017337348188092679, 6.103036322765087)     # dwave.samplers.sa.sampler._default_ising_beta_range on K2000
APC_LAMBDA = 88.9057908343783                               # largest eigenvalue of -J (= w) for WK2000_1 (selftest re-derives)
STATICA_T_INIT = (30.0, 40.0, 50.0)
STATICA_Q, STATICA_T_FIN = 4.0, 5.0
TEC_JV, TEC_T0, TEC_T1 = 30.0, 100.0, 0.1
APC_RQ, APC_QLIM, APC_T0, APC_T1 = 0.45, 0.0, 10.0, 0.1
REAIM_KSET, REAIM_T1 = (128, 256, 512, 1024), 0.1
ASB_DT, ASB_M, ASB_XI = 0.9, 2, 1.0
SB_DTS, SB_XI = (0.25, 0.5, 0.75, 1.0, 1.25), 1.0


def energy_from_field(s, h):
    return -0.5 * np.einsum('bi,bi->b', s, h)


# ------------------------------------------------------------------ TEC as published: Glauber SA with temporal coupling
@nb.njit(cache=True)
def _tec_glauber_kernel(J, s, Tsched, Jv, seeds, E):
    """Sequential heat-bath sweeps (index order). Cycle t uses the temporal field Jv*prev_i, prev = configuration at
    the start of cycle t (end of cycle t-1; the start state for t = 0). Flip probability 1/(1 + exp(2 beta s_i
    (h_i + Jv prev_i))) (TEC Eq. 4). One uniform per update attempt. E[t] = H after t cycles."""
    B, N = s.shape
    S = Tsched.shape[0]
    for b in range(B):
        np.random.seed(seeds[b])
        sb = s[b].copy()
        h = np.zeros(N)
        for i in range(N):
            acc = 0.0
            for j in range(N):
                acc += J[i, j] * sb[j]
            h[i] = acc
        e = 0.0
        for i in range(N):
            e += sb[i] * h[i]
        E[0, b] = -0.5 * e
        prev = sb.copy()
        for t in range(S):
            beta = 1.0 / Tsched[t]
            for i in range(N):
                prev[i] = sb[i]
            for i in range(N):
                x = 2.0 * beta * sb[i] * (h[i] + Jv * prev[i])
                u = np.random.random()
                if x > 700.0:
                    pf = 0.0
                elif x < -700.0:
                    pf = 1.0
                else:
                    pf = 1.0 / (1.0 + math.exp(x))
                if u < pf:
                    sb[i] = -sb[i]
                    d = 2.0 * sb[i]
                    for j in range(N):
                        h[j] += d * J[i, j]
            e = 0.0
            for i in range(N):
                e += sb[i] * h[i]
            E[t + 1, b] = -0.5 * e
        s[b] = sb


def geom_sched(T0, T1, S):
    return T0 * (T1 / T0) ** (np.arange(S) / max(1, S - 1))


def sa_neal_chunk(J, s0, seeds, S):
    """Neal defaults: beta geometric from NEAL_BETA[0] to NEAL_BETA[1] (np.geomspace), i.e. T geometric 1/beta."""
    Tsched = 1.0 / np.geomspace(NEAL_BETA[0], NEAL_BETA[1], num=S)
    s = s0.astype(np.float64).copy(); E = np.zeros((S + 1, s.shape[0]))
    M._sa_kernel(J.astype(np.float64), s, Tsched.astype(np.float64), np.asarray(seeds, np.int64), E)
    return s.astype(np.float32), E


def tec_pub_chunk(J, s0, seeds, S):
    Tsched = geom_sched(TEC_T0, TEC_T1, S)
    s = s0.astype(np.float64).copy(); E = np.zeros((S + 1, s.shape[0]))
    _tec_glauber_kernel(J.astype(np.float64), s, Tsched.astype(np.float64), float(TEC_JV), np.asarray(seeds, np.int64), E)
    return s.astype(np.float32), E


# ------------------------------------------------------------------ APC-SCA as published (Algorithm 2)
def apc_pub_traj(J, s0, S, rng, lam=APC_LAMBDA, r_q=APC_RQ, q_limit=APC_QLIM, T0=APC_T0, T1=APC_T1):
    T = geom_sched(T0, T1, S)                       # Eq. (3): T(s) = T_init r_T^(s-1), r_T = (T_final/T_init)^(1/(S-1))
    s = s0.astype(np.float32).copy(); h = s @ J
    B, N = s.shape
    E = np.empty((S + 1, B)); E[0] = energy_from_field(s, h)
    q0 = float(lam) / 2
    q = np.full((B, N), q0, np.float64)
    for t, Tt in enumerate(T):
        z = (s * h).astype(np.float64) + q                                 # h s is an exact integer in float32
        pflip = 1.0 / (1.0 + np.exp(np.clip(z / Tt, -700.0, 700.0)))     # sigmoid(-(h s + q_i)/T)
        flip = pflip > rng.random(s.shape, dtype=np.float32)              # Algorithm 2 line 8: if p_i > rand
        q = np.where(flip, q0, np.maximum(q * r_q, q_limit))
        d = np.where(flip, -2 * s, 0).astype(np.float32)
        s += d; h += d @ J
        E[t + 1] = energy_from_field(s, h)
    return s, E


# ------------------------------------------------------------------ configurations
def cfg_of(method, S, k=0):
    """Configuration k of a baseline at budget S (k indexes the paper's own candidate set where one exists)."""
    if method == 'SA (Neal)':
        return dict(family='SA_neal', S=S, beta_range=list(NEAL_BETA))
    if method == 'SCA (STATICA)':
        return dict(family='plain', q=STATICA_Q, T0=STATICA_T_INIT[k], S=S)
    if method == 'TEC on STATICA':
        return dict(family='tec', jv=TEC_JV, q=STATICA_Q, T0=STATICA_T_INIT[k], S=S)
    if method == 'TEC (published)':
        return dict(family='TEC_glauber', jv=TEC_JV, T0=TEC_T0, T1=TEC_T1, S=S)
    if method == 'APC-SCA (published)':
        return dict(family='APC_pub', lam=APC_LAMBDA, r_q=APC_RQ, q_limit=APC_QLIM, T0=APC_T0, T1=APC_T1, S=S)
    if method == 'ReAIM ASA':
        return dict(family='ReAIM', kset=REAIM_KSET, T1=REAIM_T1, S=S)
    if method == 'aSB':
        return dict(family='aSB', dt=ASB_DT, xi=ASB_XI, M=ASB_M, S=S)
    if method in ('bSB', 'dSB'):
        return dict(family=method, dt=SB_DTS[k], xi=SB_XI, S=S)
    raise ValueError(method)


N_CFG = {'SA (Neal)': 1, 'SCA (STATICA)': 3, 'TEC on STATICA': 3, 'TEC (published)': 1, 'APC-SCA (published)': 1,
         'ReAIM ASA': 1, 'aSB': 1, 'bSB': 5, 'dSB': 5}
CODE = {'SA (Neal)': 1, 'SCA (STATICA)': 2, 'TEC (published)': 3, 'TEC on STATICA': 4, 'APC-SCA (published)': 5,
        'ReAIM ASA': 6, 'aSB': 7, 'bSB': 8, 'dSB': 9, 'TEC-T (ours)': 10, 'Onsager SCA (ours)': 11}
SEQUENTIAL = ('SA_neal', 'TEC_glauber')


def draw_sequential_inputs(seed, B, N):
    """Same RNG consumption as methods.run_method for SA: initial spins, then one numba seed per run."""
    rng = np.random.default_rng(seed)
    s0 = rng.choice(np.array([-1.0, 1.0], dtype=np.float32), size=(B, N))
    seeds = rng.integers(0, 2 ** 31 - 1, size=B)
    return s0, seeds


def run_vectorized(J, sumw, cfg, B, seed):
    """All non-sequential methods. Returns E (S+1, B), final cuts, flips (SCA family only, else None)."""
    fam = cfg['family']; S = cfg['S']
    if fam == 'APC_pub':
        rng = np.random.default_rng(seed)
        s0 = rng.choice(np.array([-1.0, 1.0], dtype=np.float32), size=(B, len(J)))
        s, E = apc_pub_traj(J, s0, S, rng, cfg['lam'], cfg['r_q'], cfg['q_limit'], cfg['T0'], cfg['T1'])
        return E, (sumw - E[-1]) / 2.0, None
    c = {k: v for k, v in cfg.items() if k != 'beta_range'}
    E, cuts = M.run_method(J, sumw, c, B, seed)
    return E, cuts, None


def run_pilot_sca(J, sumw, cfg, B, seed):
    """Pilot for the SCA-family selections: abl.run (identical arithmetic to methods.sca_traj) returns flips per run."""
    rng = np.random.default_rng(seed)
    s0 = rng.choice(np.array([-1.0, 1.0], dtype=np.float32), size=(B, len(J)))
    s, flips = abl.run(J, s0, dict(cfg), rng)
    return m.cut(J, sumw, s), flips


def statica_t_ms(mean_flips, S):
    """STATICA time model fitted to the published Fig. 22 times (research/statica_reproduction_20261003): cycles at 300 MHz."""
    return (0.434 * mean_flips - 0.98 * S + 1989) / 300e3


def r99(p):
    return math.inf if p <= 0 else (1.0 if p >= 1 else math.log(0.01) / math.log1p(-p))


def mcs99(S, p):
    return math.inf if p <= 0 else (S if p >= 1 else S * math.log(0.01) / math.log1p(-p))


def wilson(k, n, z=1.959963984540054):
    p = k / n; d = 1 + z * z / n
    c = (p + z * z / (2 * n)) / d; h = z * math.sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / d
    return max(0.0, c - h), min(1.0, c + h)


def summarize_final(E, cuts, sumw, S):
    """The frozen sweep.py final record, plus best-visited statistics (secondary estimator)."""
    E = np.asarray(E, np.float64); cuts = np.asarray(cuts, np.float64); B = len(cuts)
    valid = bool(np.isfinite(cuts).all())
    k = int((cuts >= 33000).sum()); p = k / B
    best = (sumw - np.nanmin(E, axis=0)) / 2.0
    kb = int((best >= 33000).sum())
    return dict(runs=B, valid=valid, mean_cut=float(cuts.mean()), sd_cut=float(cuts.std(ddof=1)), max_cut=float(cuts.max()),
                k33000=k, p33000=p, p33000_wilson95=wilson(k, B), p33200=float((cuts >= 33200).mean()),
                p33337=float((cuts >= 33337).mean()), mcs99=mcs99(S, p),
                H_mean=E.mean(1).tolist(), H_p10=np.percentile(E, 10, axis=1).tolist(), H_p90=np.percentile(E, 90, axis=1).tolist(),
                cuts=cuts.tolist(), best_visited=dict(mean_cut=float(best.mean()), k33000=kb, p33000=kb / B,
                                                      p33000_wilson95=wilson(kb, B), mcs99=mcs99(S, kb / B)))


# ------------------------------------------------------------------ self-test
def _literal_tec(J, s0, Tsched, Jv, u):
    """Scalar transcription of TEC Eq. (4) with supplied uniforms u[t, i] (for checking the numba kernel)."""
    s = s0.copy(); N = len(s)
    for t in range(len(Tsched)):
        prev = s.copy(); beta = 1.0 / Tsched[t]
        for i in range(N):
            hs = sum(J[i, j] * s[j] for j in range(N))
            pf = 1.0 / (1.0 + math.exp(2 * beta * s[i] * (hs + Jv * prev[i])))
            if u[t, i] < pf:
                s[i] = -s[i]
    return s


def _literal_apc(J, s0, S, u, lam, r_q, q_limit, T0, T1):
    """Scalar transcription of APC-SCA Algorithm 2 with supplied uniforms u[s, i]."""
    N = len(s0); s = s0.astype(np.float64).copy(); q = np.full(N, lam / 2)
    T = geom_sched(T0, T1, S)
    for t in range(S):
        tau = s.copy(); qn = q.copy()
        for i in range(N):
            h = sum(J[i, j] * s[j] for j in range(N))
            p = 1.0 / (1.0 + math.exp((h * s[i] + q[i]) / T[t]))
            if p > u[t, i]:
                tau[i] = -s[i]; qn[i] = lam / 2
            else:
                qn[i] = max(q[i] * r_q, q_limit)
        s, q = tau, qn
    return s


def selftest():
    import hashlib
    # SHA-256 of the published files (local paths were made relative for publication)
    exp = {'research/ablation_20261005/abl.py': 'cf06c2b851f5605417ad3bd28d703a84d7ab62ac38c4eadd42c1f667470f050c',
           'research/algorithm_compare_20261007/methods.py': 'd745bdb7119e042b31c04f4c26d1217158dceef7ab9e76c83b4e4131a0965ee4',
           'research/reaim_reproduction_20261003/pilot.py': '68a436074ac8d766bf15789673e32a523b92061173f732ae4a0edda57f52093b',
           'research/theory_ideas_20261003/sweeps_abc.py': '2aa701f0c2e5c3315629577e75d2d4209f1d6f89d7636cc384515211a4cd62df',
           'research/statica_reproduction_20261003/statica_repro.py': 'f3d9f101697c8521b28c472f849fd54b2620fa4df68c47bb2b637ac0c8e33bad'}
    for f, h in exp.items():
        assert hashlib.sha256((ROOT / f).read_bytes()).hexdigest() == h, f
    J, sumw = m.load()
    # APC lambda and Neal beta range
    lam = float(np.linalg.eigvalsh(-J.astype(np.float64))[-1])
    assert abs(lam - APC_LAMBDA) < 1e-6, lam
    assert abs(NEAL_BETA[0] - math.log(2) / (2 * 1999)) < 1e-15 and abs(NEAL_BETA[1] - math.log(2000 / 0.01) / 2) < 1e-12
    # TEC kernel vs literal scalar transcription on a small random graph, shared uniforms (via a seeded generator)
    rs = np.random.default_rng(5)
    for trial in range(4):
        n = 10; W = np.triu(rs.choice([-1.0, 1.0], size=(n, n)), 1); Js = -(W + W.T)
        s0 = rs.choice([-1.0, 1.0], size=n); Ts = geom_sched(8.0, 0.5, 6)
        # numba kernel with a known seed: reproduce its uniforms by re-seeding numba's legacy generator in Python
        sk = s0[None, :].copy(); Ek = np.zeros((7, 1))
        _tec_glauber_kernel(Js, sk, Ts, 1.5, np.array([1234 + trial], np.int64), Ek)
        u = _numba_uniforms(1234 + trial, 6 * n).reshape(6, n)
        lit = _literal_tec(Js, s0.copy(), Ts, 1.5, u)
        assert np.array_equal(sk[0], lit), 'TEC kernel differs from the literal transcription'
    # APC vectorized vs literal Algorithm 2 (shared uniforms through a float32 generator)
    for trial in range(3):
        n = 9; W = np.triu(rs.choice([-1.0, 1.0], size=(n, n)), 1); Js = (-(W + W.T)).astype(np.float32)
        s0 = rs.choice(np.array([-1.0, 1.0], dtype=np.float32), size=(1, n)); lam_s = float(np.linalg.eigvalsh(-Js.astype(np.float64))[-1])
        g = np.random.default_rng(77 + trial); sv, Ev = apc_pub_traj(Js, s0, 7, g, lam_s, 0.45, 0.0, 10.0, 0.1)
        g = np.random.default_rng(77 + trial); u = np.stack([g.random((1, n), dtype=np.float32)[0] for _ in range(7)])
        lit = _literal_apc(Js.astype(np.float64), s0[0].astype(np.float64), 7, u.astype(np.float64), lam_s, 0.45, 0.0, 10.0, 0.1)
        assert np.array_equal(sv[0].astype(np.float64), lit), 'APC differs from the literal transcription'
    # Neal path = methods.sa_traj with T0 = 1/beta_hot, T1 = 1/beta_cold (geometric in T == geometric in beta)
    s0, seeds = draw_sequential_inputs(np.random.SeedSequence([1, 2, 3]), 3, len(J))
    a = sa_neal_chunk(J, s0, seeds, 9)[1]
    b = M.sa_traj(J, s0, dict(T0=1 / NEAL_BETA[0], T1=1 / NEAL_BETA[1], S=9), seeds)[1]
    assert np.allclose(a, b) and np.array_equal(a[-1], b[-1]), 'Neal schedule differs from sa_traj'
    # chunking does not change sequential results
    s0, seeds = draw_sequential_inputs(np.random.SeedSequence([4, 5, 6]), 4, len(J))
    full = tec_pub_chunk(J, s0, seeds, 5)[1]
    part = np.concatenate([tec_pub_chunk(J, s0[:2], seeds[:2], 5)[1], tec_pub_chunk(J, s0[2:], seeds[2:], 5)[1]], axis=1)
    assert np.array_equal(full, part)
    print('pub_methods selftest OK  lambda =', lam)


def _numba_uniforms(seed, n):
    @nb.njit
    def draw(seed, n):
        np.random.seed(seed)
        out = np.empty(n)
        for i in range(n):
            out[i] = np.random.random()
        return out
    return draw(seed, n)


if __name__ == '__main__':
    selftest()
