"""Independent re-implementations for headline numbers of the published-settings Figure 2 (verification only; fresh seeds
SeedSequence([20261008, 621, case]); not part of any table):
  bsb_250   bSB written from Goto et al. 2021 Eqs. (15)-(16) in float64 numpy (independent of methods.sb_traj):
            a0 = 1, a(t_k) = k/S, c0 = 0.5/(<J> sqrt N), Delta t = 1.25, walls, x, y ~ U(-0.1, 0.1); S = 250, 512 runs
            -> P(cut >= 33,000) (PROTOCOL_PUB: 254/256 = 0.992, MCS99 237)
  statica_4000  STATICA via research/statica_reproduction_20261003/statica_repro.sca (an independent implementation of
            Algorithm 1 with Eq. 7), q 4, T 40 -> 5, S 4000, 256 runs (PROTOCOL_PUB T_init 40 final: p 0.992, MCS99 3,796)
  sca_2000   the same at S 2000 (PROTOCOL_PUB: p 0.902)
Scoring of every final state from the raw edge list (statica_repro.load_graph edges), independent of the energy code."""
import json
import math
import sys
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
from common import M  # noqa: E402

HERE = Path(__file__).resolve().parent
sr = sys.modules['statica_repro']


def cut_edges(s, e):
    i, j, w = e
    return ((s[:, i] != s[:, j]) * w[None, :]).sum(1)


def bsb(J, B, S, dt, rng):
    N = J.shape[0]
    Jd = J.astype(np.float64)
    sig = math.sqrt((Jd ** 2).sum() / (N * (N - 1)))
    c0 = 0.5 / (sig * math.sqrt(N))
    x = rng.uniform(-0.1, 0.1, size=(B, N)); y = rng.uniform(-0.1, 0.1, size=(B, N))
    for k in range(S):
        a = k / S
        y = y + (-(1.0 - a) * x + c0 * (x @ Jd)) * dt
        x = x + y * dt
        out = np.abs(x) > 1
        x = np.where(out, np.sign(x), x); y = np.where(out, 0.0, y)
    return np.where(x >= 0, 1.0, -1.0)


def job(case):
    Jf, edges = sr.load_graph(); J = Jf.astype(np.float32)
    e = (edges[0], edges[1], edges[2])
    rng = np.random.default_rng(np.random.SeedSequence([20261008, 621, {'bsb_250': 0, 'statica_4000': 1, 'sca_2000': 2}[case]]))
    if case == 'bsb_250':
        s = bsb(J, 512, 250, 1.25, rng)
    else:
        S = 4000 if case == 'statica_4000' else 2000
        B = 256
        s0 = rng.choice(np.array([-1.0, 1.0], dtype=np.float32), size=(B, len(J)))
        s, _ = sr.sca(J, S, 40.0, 5.0, 4.0, s0, lambda t: rng.random((B, len(J)), dtype=np.float32))
    c = cut_edges(np.asarray(s), e)
    k = int((c >= 33000).sum())
    return dict(case=case, runs=len(c), k=k, p=k / len(c), mean_cut=float(c.mean()))


def main():
    with ProcessPoolExecutor(max_workers=3) as ex:
        out = list(ex.map(job, ['bsb_250', 'statica_4000', 'sca_2000']))
    for r in out:
        print(r, flush=True)
    (HERE / 'results_validation' / 'indep_check.json').write_text(json.dumps(out, indent=1) + '\n')


if __name__ == '__main__':
    main()
