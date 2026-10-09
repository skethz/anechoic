"""others.py against research/algorithm_compare_20261007/methods.py (unedited), on synthetic graphs only.
SA, ReAIM: identical final spins expected (sparse graph). dSB: identical expected on a complete +-1 graph (sigma_J = 1).
bSB, aSB: real-valued J.x summed in a different order -> compare mean final cut statistically (complete +-1 graph,
where methods.py's normalisation equals ours)."""
import json
import math
import sys

import numpy as np

import engine as E
import others as O
import synth
import methods as M   # noqa: E402  (path set up by engine)


def complete_pm(n, seed):
    rng = np.random.default_rng(seed)
    iu = np.triu_indices(n, 1)
    return E.Graph(f'syn_complete_{n}_{seed}', n, (iu[0], iu[1], rng.choice(np.array([-1, 1]), size=len(iu[0]))))


def main():
    out = {}; ok = True
    g = synth.random_graph(300, 2700, True, 11); J = g.dense()
    # SA
    cfg = dict(family='SA', T0=30 * g.alpha, T1=1.0 * g.alpha, S=60)
    rng = np.random.default_rng(5); s0 = rng.choice(np.array([-1.0, 1.0], dtype=np.float32), size=(8, g.N))
    sref, _ = M.sa_traj(J, s0, cfg, rng.integers(0, 2 ** 31 - 1, size=8))
    snew = O.sa(g, 8, cfg, np.random.default_rng(5))
    out['SA_identical'] = bool(np.array_equal(sref, snew)); ok &= out['SA_identical']
    # ReAIM
    cfg = dict(family='ReAIM', kset=(4, 8, 16, 32), T1=0.1, S=300)
    sref, _ = M.reaim_traj(J, 6, cfg, np.random.default_rng(9))
    snew = O.reaim(g, 6, cfg, np.random.default_rng(9))
    out['ReAIM_identical'] = bool(np.array_equal(sref, snew)); ok &= out['ReAIM_identical']
    # SB family on a complete +-1 graph (sigma_J = 1)
    gc = complete_pm(200, 12); Jc = gc.dense(); W = float(gc.W)
    assert abs(O.sigma_J(gc) - 1.0) < 1e-12
    cfg = dict(family='dSB', dt=1.0, xi=1.0, S=60)
    sref, _ = M.sb_traj(Jc, 6, cfg, np.random.default_rng(3)); snew = O.sb(gc, 6, cfg, np.random.default_rng(3))
    out['dSB_identical_S60'] = bool(np.array_equal(sref, snew))
    for fam, cfg in (('bSB', dict(family='bSB', dt=1.0, xi=1.0, S=400)), ('dSB', dict(family='dSB', dt=1.0, xi=1.0, S=400)),
                     ('aSB', dict(family='aSB', dt=0.9, xi=1.0, S=400))):
        B = 128
        if fam == 'aSB':
            sref, _, fin = M.asb_traj(Jc, B, cfg, np.random.default_rng(21)); snew, fin2 = O.asb(gc, B, cfg, np.random.default_rng(22))
        else:
            sref, _ = M.sb_traj(Jc, B, cfg, np.random.default_rng(21)); snew = O.sb(gc, B, cfg, np.random.default_rng(22))
        cr, cn = gc.cut(sref), gc.cut(snew)
        se = math.sqrt(cr.var(ddof=1) / B + cn.var(ddof=1) / B)
        zsc = (cn.mean() - cr.mean()) / se
        out[f'{fam}_mean_cut_ref'] = float(cr.mean()); out[f'{fam}_mean_cut_new'] = float(cn.mean()); out[f'{fam}_z'] = float(zsc)
        ok &= abs(zsc) < 3.5
        print(fam, 'ref', cr.mean(), 'new', cn.mean(), 'z', round(zsc, 2), flush=True)
    out['all_pass'] = bool(ok)
    print(json.dumps(out, indent=1))
    (E.ROOT / 'verify_others.json').write_text(json.dumps(out, indent=1) + '\n')
    return 0 if ok else 1


if __name__ == '__main__':
    sys.exit(main())
