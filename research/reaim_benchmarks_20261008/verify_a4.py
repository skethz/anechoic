"""Checks for Addendum 4 kernels before any A4 run (writes verify_a4_<host>.json):
 1. solvers_a4.tec_seq equals a pure-Python transcription of the sequential Glauber TEC update (legacy MT stream per
    run, same draws) in final spins and flips, for J_v in {0, 30, -6}, on small Max-Cut, GPP (uniform coupling) and TSP
    (bias) instances; the energy trace's last value equals the energy recomputed from the final spins.
 2. solvers_a4.sca_sig_q equals a numpy transcription of Okonogi 2023 Algorithm 1 (synchronous exact sigmoid,
    exponential T and q, argmin output)."""
import json
import math
import os
import platform
import sys
from pathlib import Path

sys.dont_write_bytecode = True
HERE = Path(__file__).resolve().parent
os.environ.setdefault('NUMBA_CACHE_DIR', str(HERE / '.numba_cache'))
sys.path.insert(0, str(HERE))
import numpy as np  # noqa: E402

import problems as PB  # noqa: E402
import solvers as SV  # noqa: E402
import solvers_a4 as SV4  # noqa: E402


def tec_ref(P, B, cfg, rng):
    s0 = rng.choice(np.array([-1.0, 1.0], dtype=np.float32), size=(B, P.N)).astype(np.float64)
    seeds = rng.integers(0, 2 ** 31 - 1, size=B)
    J = P.dense(); b = P.b.astype(np.float64)
    T = SV.sched(cfg['T0'], cfg['S'], cfg['tfin'])
    out = s0.copy(); flips = np.zeros(B)
    for r in range(B):
        np.random.seed(int(seeds[r]))
        s = s0[r].copy(); prev = s.copy()
        for Tt in T:
            for i in range(P.N):
                hs = J[i] @ s + b[i]
                a = 2.0 / Tt * s[i] * (hs + cfg['jv'] * prev[i])
                pf = 0.0 if a > 700 else (1.0 if a < -700 else 1.0 / (1.0 + math.exp(a)))
                if np.random.random() < pf:
                    s[i] = -s[i]; flips[r] += 1
            prev = s.copy()
        out[r] = s
    return out.astype(np.float32), flips


def sca_q_ref(P, B, cfg, rng):
    s = rng.choice(np.array([-1.0, 1.0], dtype=np.float32), size=(B, P.N)).astype(np.float64)
    J = P.dense(); b = P.b.astype(np.float64)
    S = cfg['S']; T = SV.sched(cfg['T0'], S, cfg['tfin'])
    qs = cfg['q_init'] * (cfg['q_final'] / cfg['q_init']) ** (np.arange(S) / max(1, S - 1))
    E = lambda x: -0.5 * np.einsum('bi,ij,bj->b', x, J, x) - x @ b  # noqa: E731
    best = s.copy(); bE = E(s)
    for t, Tt in enumerate(T):
        u = rng.random(s.shape)
        z = s * (s @ J + b[None, :]) + qs[t]
        pf = 1.0 / (1.0 + np.exp(np.clip(z / Tt, -700, 700)))
        s = np.where(u < pf, -s, s)
        e = E(s); imp = e < bE
        best[imp] = s[imp]; bE[imp] = e[imp]
    return best.astype(np.float32)


def main():
    out = dict(checks=[]); ok = True
    rng0 = np.random.default_rng(13)
    n = 30; Ed = [(i, j) for i in range(n) for j in range(i + 1, n) if rng0.random() < 0.25]
    ei = np.array([e[0] for e in Ed]); ej = np.array([e[1] for e in Ed]); w = rng0.choice(np.array([-1, 1]), size=len(Ed))
    mc = PB.Problem('v_mcp30', 'mcp', n, np.r_[ei, ej], np.r_[ej, ei], -np.r_[w, w].astype(float), 0.0, np.zeros(n), 1.0,
                    dict(ei=ei, ej=ej, w=w, ref=1, target=1))
    gp = PB.gpp_from_edges('v_gpp30', n, ei, ej, np.ones(len(Ed), np.int64), 4, ref=1, target=1)
    pts = rng0.random((4, 2)) * 100
    W = np.rint(np.hypot(*(pts[:, None, :] - pts[None, :, :]).transpose(2, 0, 1))).astype(np.int64); np.fill_diagonal(W, 0)
    W[(W == 0) & ~np.eye(4, dtype=bool)] = 1
    tp = PB.tsp_from_matrix('v_tsp4', W, float(W.max()), ref=1, target=1)
    for P in (mc, gp, tp):
        a = P.sigma_T
        for jv in (0.0, 30.0, -6.0):
            cfg = dict(family='tec_seq', jv=jv * a / 10, T0=10.0 * a, tfin=0.1 * a, S=40)
            s1, f1, Etr = SV4.tec_seq(P, 3, cfg, np.random.default_rng(np.random.SeedSequence([3, 3])), trace=True)
            s2, f2 = tec_ref(P, 3, cfg, np.random.default_rng(np.random.SeedSequence([3, 3])))
            Efin = P.energy(s1.astype(np.float64))
            same = bool(np.array_equal(s1, s2) and np.array_equal(f1, f2) and np.allclose(Etr[:, -1], Efin))
            ok &= same
            out['checks'].append(dict(check='tec_seq == transcription', problem=P.name, jv=jv, identical=same, flips=f1.tolist()))
            print('tec_seq', P.name, jv, same, flush=True)
        cfg = dict(q_init=1.5 * a, q_final=0.2 * a, T0=10.0 * a, tfin=0.1 * a, S=60)
        b1, _, _ = SV4.sca_sig_q(P, 3, cfg, np.random.default_rng(np.random.SeedSequence([4])))
        b2 = sca_q_ref(P, 3, cfg, np.random.default_rng(np.random.SeedSequence([4])))
        same = bool(np.array_equal(b1, b2)); ok &= same
        out['checks'].append(dict(check='sca_sig_q == transcription', problem=P.name, identical=same))
        print('sca_sig_q', P.name, same, flush=True)
    out['all_ok'] = bool(ok)
    out['platform'] = dict(node=platform.node(), machine=platform.machine())
    (HERE / f"verify_a4_{platform.node().split('.')[0]}.json").write_text(json.dumps(out, indent=1) + '\n')
    print('ALL OK' if ok else 'FAIL')
    return 0 if ok else 1


if __name__ == '__main__':
    sys.exit(main())
