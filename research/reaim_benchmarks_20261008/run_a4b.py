"""Addendum 4.1: the SB methods on TSP with the papers' treatment of local fields (Goto 2021: ancillary spin), fresh seeds,
fresh directory results_a4/final/sb_anc/. Variants:
  bSB, dSB: a0 = 1, c0 = 0.5/(<J'> sqrt(N+1)) on the (N+1)-spin field-free problem; dt by the paper's procedure (best of
            {0.25, 0.5, 0.75, 1, 1.25} by 64-run pilots, ties to the smaller dt); seeds pilot [.., 55, p, i, m, s, k],
            final [.., 54, p, i, m, s]
  aSB_dt09: Goto 2019 K2000 values (dt = 0.9, M = 2), xi0 = 0.7/(SD(J') sqrt(N+1)); seeds [.., 56, ..]
  aSB_dt05: dt = 0.5 (OUR CHOICE, as in Addendum 4), M = 2; seeds [.., 57, ..]
Usage: python3 run_a4b.py run [W]; python3 run_a4b.py check"""
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

OUT = HERE / 'results_a4' / 'final' / 'sb_anc'
VARIANTS = {'bSB': ('bSB', 54), 'dSB': ('dSB', 54), 'aSB_dt09': ('aSB', 56), 'aSB_dt05': ('aSB', 57)}


def run_sb(Q, cfg, B, seed):
    rng = np.random.default_rng(seed)
    if cfg['family'] == 'aSB':
        s, finite = SV.asb(Q, B, cfg, rng)
        return AN.decode(s), finite
    return AN.decode(SV.sb(Q, B, cfg, rng)), None


def job(args):
    variant, inst, S = args
    method, sk = VARIANTS[variant]
    path = OUT / variant / RA.tag('tsp', inst) / f'{method}_S{S}.json'
    if path.exists():
        return None
    t0 = time.time()
    P = RA.problem('tsp', inst); Q = AN.with_ancilla(P)
    C = ST.instance_constants(Q)
    mi = R3.METHODS.index(method); si = GA.S_LIST['tsp'].index(S); pc = RA.PCODE['tsp']; ic = RA.icode('tsp', inst)
    pilot = None
    if method in ('bSB', 'dSB'):
        pilot = []
        for k, dt in enumerate(ST.SB_DT_SET):
            c = dict(family=method, dt=dt, xi=1.0, S=S)
            s, _ = run_sb(Q, c, R3.B_PILOT, np.random.SeedSequence([R3.SEED, 55, pc, ic, mi, si, k]))
            pilot.append(dict(dt=dt, mean_quality=float(P.evaluate(s)['quality'].mean())))
        sel = max(range(len(pilot)), key=lambda k: (pilot[k]['mean_quality'], -k))
        cfg = dict(family=method, dt=ST.SB_DT_SET[sel], xi=1.0, S=S, note='Goto 2021 with ancillary spin; dt by the paper procedure')
    else:
        dt = 0.9 if variant == 'aSB_dt09' else 0.5
        cfg = dict(family='aSB', dt=dt, M=2, xi=C['rms_J'] / C['sd_J'], S=S,
                   note=f'Goto 2019 with ancillary spin; dt = {dt}' + (' (OUR CHOICE)' if dt == 0.5 else ' (paper K2000 value)'))
    s, finite = run_sb(Q, cfg, R3.B_FINAL, np.random.SeedSequence([R3.SEED, sk, pc, ic, mi, si]))
    ev = P.evaluate(s)
    fin, q, feas = R3.stats(P, 'tsp', ev, None, finite, S)
    rec = dict(variant=variant, problem='tsp', instance=P.name, method=method, S=S, cfg=cfg, dt_pilot=pilot, final=fin,
               quality=q.tolist(), values=ev['value'].tolist(), feasible_runs=feas.tolist(),
               finite=None if finite is None else np.asarray(finite).tolist(), sec=time.time() - t0)
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_suffix('.tmp'); tmp.write_text(json.dumps(rec, default=list) + '\n'); tmp.rename(path)
    return (f"{variant} {P.name} S={S} dt={cfg['dt']} q={fin['mean_quality']:.4f} feas={fin['p_feasible']:.3f} "
            f"p={fin['p_target']:.3f} invalid={fin['n_invalid']} ({rec['sec']:.0f}s)")


def check():
    """The ancilla problem reproduces the energy: E'(s, s0 = +1) == E(s) and E'(-s, -1) == E(s) on random states."""
    ok = True
    for inst in GA.TSP_INST:
        P = RA.problem('tsp', inst); Q = AN.with_ancilla(P)
        rng = np.random.default_rng(1)
        s = rng.choice(np.array([-1.0, 1.0]), size=(8, P.N))
        sp = np.concatenate([s, np.ones((8, 1))], 1); sm = -sp
        same = np.allclose(Q.energy(sp), P.energy(s)) and np.allclose(Q.energy(sm), P.energy(s)) and \
            np.array_equal(AN.decode(sm), s.astype(np.float32))
        ok &= bool(same)
        print(inst, 'ancilla energy/decoding identical:', bool(same), flush=True)
    print('ALL OK' if ok else 'FAIL')
    return ok


def main():
    W = int(sys.argv[2]) if len(sys.argv) > 2 else 64
    J = [(v, inst, S) for v in VARIANTS for inst in GA.TSP_INST for S in GA.S_LIST['tsp']]
    J.sort(key=lambda j: -(RA.problem('tsp', j[1]).N * j[2] * (5 if j[0] in ('bSB', 'dSB') else 1)))
    RA.run_pool(J, job, W, 'A4.1 SB on TSP with ancillary spin')


if __name__ == '__main__':
    if sys.argv[1] == 'check':
        sys.exit(0 if check() else 1)
    main()
