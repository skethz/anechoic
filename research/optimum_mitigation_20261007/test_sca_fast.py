"""Validation of sca_fast against abl.run (development check; not part of the pre-registered measurements).
1. Exact replay: identical uniforms injected into abl.run and into sca_fast.run_steps -> identical final spins and flips.
2. Speed: trials per second at S = 8000 (Onsager-online) with the thread count in NUMBA_NUM_THREADS.
Usage: python3 test_sca_fast.py [exact|speed|stat]"""
import json
import sys
import time

import numpy as np

import sca_fast as F


class FakeRng:
    def __init__(self, U):
        self.U = U; self.k = 0

    def random(self, shape, dtype=np.float32):
        u = self.U[self.k]; self.k += 1
        assert u.shape == tuple(shape)
        return u


def exact():
    J, sumw, J8 = F.load(); J2 = (2 * J8).astype(np.int8)
    cfgs = [dict(family='plain', q=8.0, T0=30.0, S=300),
            dict(family='onsager', q=8.0, lam=1.05, ramp=True, T0=12.0, S=300),
            dict(family='onsager', q=6.0, lam=0.9, ramp=False, T0=12.0, S=200),
            dict(family='tecT', q=8.0, kappa=1.75, ramp=True, T0=15.0, S=280)]
    out = []
    for i, cfg in enumerate(cfgs):
        B = 6
        g = np.random.default_rng([777, i])
        s0 = np.where(g.random((B, 2000)) < 0.5, -1, 1).astype(np.int8)
        U = g.random((cfg['S'], B, 2000), dtype=np.float32)
        s_ref, fl_ref = F.abl.run(J, s0.astype(np.float32), cfg, FakeRng(U))
        st = F.State(J, None, s0, np.arange(B, dtype=np.uint64))
        F.advance(J2, st, F.tables(cfg), 0, cfg['S'], U=U)
        same_s = bool(np.array_equal(st.s.astype(np.float32), s_ref))
        same_f = bool(np.array_equal(st.flips, fl_ref.astype(np.int64)))
        hf = st.s.astype(np.float32) @ J
        same_h = bool(np.array_equal(hf.astype(np.int64), st.h.astype(np.int64)))
        c_ref = F.m.cut(J, sumw, s_ref); c_new = F.cuts(st, sumw)
        r = dict(cfg=cfg, identical_spins=same_s, identical_flips=same_f, field_consistent=same_h,
                 cut_ref=c_ref.tolist(), cut_new=c_new.tolist(), flips=fl_ref.tolist())
        out.append(r); print(json.dumps(r), flush=True)
    ok = all(r['identical_spins'] and r['identical_flips'] and r['field_consistent'] for r in out)
    print('EXACT', 'PASS' if ok else 'FAIL')
    return out


def speed():
    J, sumw, J8 = F.load(); J2 = (2 * J8).astype(np.int8)
    cfg = dict(family='onsager', q=8.0, lam=1.05, ramp=True, T0=12.0, S=8000)
    F.run_cfg(J, J2, sumw, dict(cfg, S=50), 12, 1)  # compile
    for B in (24, 96):
        t = time.time(); r = F.run_cfg(J, J2, sumw, cfg, B, 2); dt = time.time() - t
        print(f'S=8000 B={B}: {dt:.1f} s, {B / dt:.2f} trials/s, mean flips {r["flips"].mean():.0f}, '
              f'mean cut {r["cut"].mean():.1f}, max {r["cut"].max():.0f}', flush=True)


def stat():
    """abl.run vs sca_fast at the probe configuration (S = 2000), 256 trials each, independent streams."""
    J, sumw, J8 = F.load(); J2 = (2 * J8).astype(np.int8)
    cfg = dict(family='onsager', q=8.0, lam=1.05, ramp=True, T0=12.0, S=2000)
    g = np.random.default_rng(4242)
    s0 = g.choice(np.array([-1.0, 1.0], dtype=np.float32), size=(256, 2000))
    t = time.time(); s, fl = F.abl.run(J, s0, cfg, g); t_ref = time.time() - t
    c_ref = F.m.cut(J, sumw, s)
    t = time.time(); r = F.run_cfg(J, J2, sumw, cfg, 256, 4243); t_new = time.time() - t
    for name, c, f, tt in (('abl.run', c_ref, fl, t_ref), ('sca_fast', r['cut'], r['flips'], t_new)):
        print(f'{name}: mean cut {c.mean():.2f} (sd {c.std():.1f}), mean flips {f.mean():.0f} (sd {f.std():.0f}), '
              f'p33000 {np.mean(c >= 33000):.3f}, p33200 {np.mean(c >= 33200):.3f}, {tt:.1f} s', flush=True)
    se = np.sqrt(c_ref.var() / 256 + r['cut'].var() / 256)
    print(f'difference of mean cut: {r["cut"].mean() - c_ref.mean():.2f} (SE {se:.2f})')


if __name__ == '__main__':
    {'exact': exact, 'speed': speed, 'stat': stat}[sys.argv[1] if len(sys.argv) > 1 else 'exact']()
