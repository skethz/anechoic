"""PROTOCOL_PUB_A2.md: ReAIM's own output (x_best, updated after every run phase, Algorithm 3) for the PROTOCOL_PUB ReAIM
finals, by re-simulating them with identical seeds. Asserts the final cuts equal the stored PROTOCOL_PUB finals.
Candidate states: the state after each run phase (t = 128, 256, ...) and the state at t = S (the last, possibly empty run
phase), plus the initial state (H_best is initialised from it in pilot.asa). Usage: python3 reaim_xbest.py [workers]"""
import json
import sys
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
import pub_methods as P  # noqa: E402
from common import M  # noqa: E402

HERE = Path(__file__).resolve().parent
S_LIST = (250, 500, 1000, 2000, 4000)
B = 256


def ends(S, it_trial=32, it_run=96):
    t, out = 0, [0]
    while t < S:
        t += min(it_trial, S - t)
        t += min(it_run, S - t)
        out.append(t)
    return sorted(set(out))


def job(si):
    J, sumw = M.m.load(); S = S_LIST[si]
    cfg = P.cfg_of('ReAIM ASA', S)
    E, cuts, _ = P.run_vectorized(J, sumw, cfg, B, np.random.SeedSequence([20261008, 302, P.CODE['ReAIM ASA'], si, 0]))
    idx = ends(S)
    best = (sumw - E[idx].min(axis=0)) / 2
    return S, cuts, best, idx


def main(workers):
    out = {}
    with ProcessPoolExecutor(max_workers=workers) as ex:
        for S, cuts, best, idx in ex.map(job, range(len(S_LIST))):
            stored = json.loads((HERE / 'results_pub' / 'raw' / f'ReAIM_ASA_S{S}_k0_final.json').read_text())['cuts']
            assert np.array_equal(np.asarray(stored), cuts), 'final cuts differ from the PROTOCOL_PUB final'
            k = int((best >= 33000).sum())
            out[S] = dict(runs=B, k33000=k, p33000=k / B, p33000_wilson95=P.wilson(k, B), mcs99=P.mcs99(S, k / B),
                          mean_cut=float(best.mean()), candidate_steps=idx, cuts_xbest=best.tolist())
            print(S, k / B, float(best.mean()), flush=True)
    (HERE / 'results_pub' / 'reaim_xbest_A2.json').write_text(json.dumps(out) + '\n')


if __name__ == '__main__':
    main(int(sys.argv[1]) if len(sys.argv) > 1 else 5)
