"""PROTOCOL_PUB_A3.md: aSB at Goto et al. 2021's K2000 aSB step (Delta t 0.5, M 2), finals for Figure 2 and a validation
against Goto 2021 Fig. 2A. Per-run divergence accounting (asb_perrun.asb_traj_flags, = methods.asb_traj arithmetic).
Usage: python3 asb_dt05.py [workers]"""
import json
import sys
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
import asb_perrun as A  # noqa: E402
import pub_methods as P  # noqa: E402
from common import M  # noqa: E402

HERE = Path(__file__).resolve().parent
S_LIST = (250, 500, 1000, 2000, 4000)
VAL = ((100, 32030.0), (1000, 32950.0), (10000, 32975.0))
JOBS = [('final', si, S) for si, S in enumerate(S_LIST)] + [('val', vi, S) for vi, (S, _) in enumerate(VAL)]


def job(j):
    np.seterr(all='ignore')
    kind, idx, S = j
    J, sumw = M.m.load()
    seed = [20261008, 304, 7, idx] if kind == 'final' else [20261008, 612, idx]
    rng = np.random.default_rng(np.random.SeedSequence(seed))
    s, E, ok = A.asb_traj_flags(J, 256, dict(family='aSB', dt=0.5, xi=1.0, M=2, S=S), rng)
    cuts = (sumw - E[-1]) / 2
    rec = dict(kind=kind, idx=idx, S=S, runs=256, diverged=int((~ok).sum()), k33000=int(((cuts >= 33000) & ok).sum()))
    good = cuts[ok]
    rec.update(mean_cut=float(good.mean()) if len(good) else None, sd_cut=float(good.std(ddof=1)) if len(good) > 1 else None)
    if kind == 'final':
        k = rec['k33000']; p = k / 256; Eg = E[:, ok]
        rec.update(p33000=p, p33000_wilson95=P.wilson(k, 256), mcs99=P.mcs99(S, p), H_mean=Eg.mean(1).tolist(),
                   H_p10=np.percentile(Eg, 10, axis=1).tolist(), H_p90=np.percentile(Eg, 90, axis=1).tolist(),
                   p33200=float(((cuts >= 33200) & ok).mean()), cuts=cuts.tolist(), ok=ok.tolist())
    else:
        ref = VAL[idx][1]; se = rec['sd_cut'] / np.sqrt(len(good))
        rec.update(paper_read=ref, pass_=bool(abs(rec['mean_cut'] - ref) <= 40 + 3 * se))
    return rec


def main(workers):
    out = {'final': {}, 'val': {}}
    with ProcessPoolExecutor(max_workers=workers) as ex:
        for r in ex.map(job, JOBS):
            out[r['kind']][str(r['S'])] = r
            print({k: v for k, v in r.items() if k not in ('H_mean', 'H_p10', 'H_p90', 'cuts', 'ok')}, flush=True)
    (HERE / 'results_pub' / 'asb_dt05_A3.json').write_text(json.dumps(out) + '\n')


if __name__ == '__main__':
    main(int(sys.argv[1]) if len(sys.argv) > 1 else 8)
