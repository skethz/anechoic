"""PROTOCOL_PUB_A1.md: aSB finals re-simulated with the same seeds and arithmetic as methods.asb_traj, plus a per-run
divergence flag (x or y non-finite at any step). Diverged runs count as failures and are excluded from the mean cut.
Checks that the non-diverged runs reproduce the PROTOCOL_PUB final cuts exactly. Usage: python3 asb_perrun.py [workers]"""
import json
import math
import sys
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
import pub_methods as P  # noqa: E402
from common import M  # noqa: E402

HERE = Path(__file__).resolve().parent
OUT = HERE / 'results_pub'
S_LIST = (250, 500, 1000, 2000, 4000)
B = 256


def asb_traj_flags(J, Bn, cfg, rng):
    """methods.asb_traj, line for line, with a per-run finiteness flag."""
    S = cfg['S']; dt = cfg['dt']; Mm = cfg.get('M', 5); N = J.shape[0]
    xi = np.float32(cfg['xi'] * 0.7 / math.sqrt(N))
    s0 = rng.choice(np.array([-1.0, 1.0], dtype=np.float32), size=(Bn, N))
    x = np.zeros((Bn, N), np.float32); y = (s0 * 0.1 * rng.random((Bn, N), dtype=np.float32)).astype(np.float32)
    E = np.empty((S + 1, Bn)); E[0] = P.energy_from_field(s0, s0 @ J)
    delta = np.float32(dt / Mm); ok = np.ones(Bn, bool); s = s0
    for t in range(1, S + 1):
        p = np.float32(t / S)
        for _ in range(Mm):
            x = x + y * delta
            y = y - ((x * x + 1 - p) * x) * delta
        y = y + (x @ J) * (xi * np.float32(dt))
        fin = np.isfinite(x).all(1) & np.isfinite(y).all(1)
        if not fin.all():
            ok &= fin; x = np.nan_to_num(x); y = np.nan_to_num(y)
        s = np.where(x >= 0, 1.0, -1.0).astype(np.float32)
        E[t] = P.energy_from_field(s, s @ J)
    return s, E, ok


def job(si):
    np.seterr(all='ignore')
    J, sumw = M.m.load(); S = S_LIST[si]
    cfg = P.cfg_of('aSB', S)
    rng = np.random.default_rng(np.random.SeedSequence([20261008, 302, P.CODE['aSB'], si, 0]))
    s, E, ok = asb_traj_flags(J, B, cfg, rng)
    cuts = (sumw - E[-1]) / 2.0
    return S, E, cuts, ok


def main(workers):
    J, sumw = M.m.load()
    out = {}
    with ProcessPoolExecutor(max_workers=workers) as ex:
        for S, E, cuts, ok in ex.map(job, range(len(S_LIST))):
            raw = json.loads(next((HERE / 'results_pub' / 'raw').glob(f'aSB_S{S}_k0_final.json')).read_text())
            # identity check: the frozen-arithmetic energies of the pre-registered final equal ours for every run
            assert np.allclose(np.array(raw['H_mean']), E.mean(1)), 'trajectory differs from the PROTOCOL_PUB final'
            k = int(((cuts >= 33000) & ok).sum()); p = k / B
            good = cuts[ok]
            Eg = E[:, ok]
            out[S] = dict(runs=B, diverged=int((~ok).sum()), k33000=k, p33000=p, p33000_wilson95=P.wilson(k, B),
                          mcs99=P.mcs99(S, p), mean_cut_nondiverged=float(good.mean()) if len(good) else None,
                          sd_cut_nondiverged=float(good.std(ddof=1)) if len(good) > 1 else None,
                          H_mean_nondiverged=Eg.mean(1).tolist() if len(good) else None,
                          H_p10_nondiverged=np.percentile(Eg, 10, axis=1).tolist() if len(good) else None,
                          H_p90_nondiverged=np.percentile(Eg, 90, axis=1).tolist() if len(good) else None)
            print(S, 'diverged', out[S]['diverged'], 'p', p, 'mean cut (non-diverged)', out[S]['mean_cut_nondiverged'], flush=True)
    (OUT / 'asb_perrun_A1.json').write_text(json.dumps(out) + '\n')


if __name__ == '__main__':
    main(int(sys.argv[1]) if len(sys.argv) > 1 else 5)
