"""Bit-equality of engine.run (sparse, fused) against abl.run (the unedited engine model) on synthetic sparse graphs.
abl.run is called with dense float32 J and abl.m.T_FIN patched to the graph's final temperature. Pass = identical final
spins and identical flip counts for every run, every family. Also checks the exact integer cut against abl.m.cut."""
import json
import sys
import time

import numpy as np

import engine as E
import synth

abl = E.abl


def check(g, cfg, tfin, B=12, seed=7):
    J = g.dense()
    rng1 = np.random.default_rng(seed); s0 = rng1.choice(np.array([-1.0, 1.0], dtype=np.float32), size=(B, g.N))
    rng2 = np.random.default_rng(seed); s0b = rng2.choice(np.array([-1.0, 1.0], dtype=np.float32), size=(B, g.N))
    abl.m.T_FIN = tfin
    t0 = time.perf_counter(); sa, fa = abl.run(J, s0, cfg, rng1); ta = time.perf_counter() - t0
    t0 = time.perf_counter(); sb, fb = E.run(g, s0b, cfg, rng2, tfin); tb = time.perf_counter() - t0
    same = bool(np.array_equal(sa, sb) and np.array_equal(fa, fb))
    cut_ref = abl.m.cut(J, float(g.W), sa); cut_new = g.cut(sb)
    return same, bool(np.array_equal(cut_ref.astype(np.int64), cut_new) and np.allclose(cut_ref, cut_new)), ta, tb, float(fa.mean())


def main():
    graphs = [synth.random_graph(300, 2700, True, 1), synth.random_graph(300, 2700, False, 2), synth.torus(16, 20, 3),
              synth.planar_like(300, True, 4), synth.planar_like(300, False, 5)]
    out = []; ok = True
    for g in graphs:
        a = g.alpha
        cfgs = [dict(family='plain', q=4.0 * a, T0=30.0 * a, S=150),
                dict(family='tec', q=6.0 * a, jv=-2.0 * a, T0=20.0 * a, S=150),
                dict(family='tec', q=6.0 * a, jv=4.0 * a, T0=20.0 * a, S=150),
                dict(family='tecT', q=8.0 * a, kappa=1.25, ramp=True, T0=20.0 * a, S=150),
                dict(family='tecT', q=6.0 * a, kappa=0.5, ramp=False, T0=12.0 * a, S=150),
                dict(family='onsager', q=8.0 * a, lam=1.05 * g.lam_factor, ramp=True, T0=12.0 * a, S=150),
                dict(family='onsager', q=6.0 * a, lam=0.5 * g.lam_factor, ramp=False, T0=15.0 * a, S=150),
                dict(family='onsager', q=6.0 * a, lam=5.0 * g.lam_factor, ramp=True, T0=15.0 * a, S=150),
                dict(family='apc', q_reset=16.0 * a, r_q=0.9, q_lim=2.0 * a, T0=20.0 * a, S=150)]
        for cfg in cfgs:
            same, cut_ok, ta, tb, fl = check(g, cfg, 5.0 * a)
            ok &= same and cut_ok
            out.append(dict(graph=g.name, cfg=cfg, identical=same, cut_identical=cut_ok, t_abl=ta, t_engine=tb, mean_flips=fl))
            print(f"{g.name:28s} {cfg['family']:8s} identical={same} cut_ok={cut_ok} abl {ta:.2f}s engine {tb:.3f}s flips {fl:.0f}",
                  flush=True)
    print('ALL IDENTICAL' if ok else 'MISMATCH')
    (E.ROOT / 'verify_engine.json').write_text(json.dumps(dict(all_identical=bool(ok), cases=out), indent=1) + '\n')
    return 0 if ok else 1


if __name__ == '__main__':
    sys.exit(main())
