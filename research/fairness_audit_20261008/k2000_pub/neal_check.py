"""PROTOCOL_PUB.md validation: the real dwave-samplers SimulatedAnnealingSampler (all defaults, num_sweeps = S) on WK2000_1,
256 reads at S = 250 and S = 1000, seed 20261008 + read index, one read per call. Compared with the kernel finals of
SA (Neal) in results_pub (pass: p and mean cut within 2.6 two-sample standard errors).
Usage: python3 neal_check.py run [workers] | compare"""
import json
import math
import sys
import time
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
from common import M  # noqa: E402

HERE = Path(__file__).resolve().parent
OUT = HERE / 'results_neal_check'
S_CHECK, READS, SEED0 = (250, 1000), 256, 20261008

_G = None


def _init():
    global _G
    import dwave.samplers as ds
    sr = sys.modules.get('statica_repro') or __import__('statica_repro')
    Jf, (i, j, w) = sr.load_graph()
    J = {(int(a), int(b)): float(c) for a, b, c in zip(i, j, w)}
    _G = (ds.SimulatedAnnealingSampler(), J, float(w.sum()), ds.__version__)


def _read(arg):
    S, r = arg
    sampler, J, W, ver = _G
    t = time.perf_counter()
    ss = sampler.sample_ising({}, J, num_reads=1, num_sweeps=S, seed=SEED0 + r)
    return S, r, (W - float(ss.first.energy)) / 2, time.perf_counter() - t, ss.info['beta_range'], ver


def run(workers):
    OUT.mkdir(exist_ok=True)
    with ProcessPoolExecutor(max_workers=workers, initializer=_init) as ex:
        res = list(ex.map(_read, [(S, r) for S in S_CHECK for r in range(READS)]))
    for S in S_CHECK:
        rows = [x for x in res if x[0] == S]
        cuts = [x[2] for x in rows]
        out = dict(S=S, reads=READS, seed0=SEED0, version=rows[0][5], beta_range=list(rows[0][4]), cuts=cuts,
                   p33000=float(np.mean(np.array(cuts) >= 33000)), mean_cut=float(np.mean(cuts)), sd_cut=float(np.std(cuts, ddof=1)),
                   mean_sec=float(np.mean([x[3] for x in rows])))
        (OUT / f'neal_S{S}.json').write_text(json.dumps(out, indent=1) + '\n')
        print(S, out['p33000'], out['mean_cut'], out['beta_range'], out['mean_sec'])


def compare():
    rows = []
    for S in S_CHECK:
        a = json.loads((OUT / f'neal_S{S}.json').read_text())
        b = json.loads((HERE / 'results_pub' / f'S{S}.json').read_text())['SA (Neal)']['final']
        n1, n2 = a['reads'], b['runs']
        p1, p2 = a['p33000'], b['p33000']; pp = (p1 * n1 + p2 * n2) / (n1 + n2)
        zp = (p1 - p2) / math.sqrt(max(1e-12, pp * (1 - pp) * (1 / n1 + 1 / n2)))
        zm = (a['mean_cut'] - b['mean_cut']) / math.sqrt(a['sd_cut'] ** 2 / n1 + b['sd_cut'] ** 2 / n2)
        rows.append(dict(S=S, neal_p=p1, kernel_p=p2, z_p=zp, neal_mean=a['mean_cut'], kernel_mean=b['mean_cut'], z_mean=zm,
                         pass_=abs(zp) <= 2.6 and abs(zm) <= 2.6))
        print(rows[-1])
    (OUT / 'compare.json').write_text(json.dumps(rows, indent=1) + '\n')


if __name__ == '__main__':
    if sys.argv[1] == 'run':
        run(int(sys.argv[2]) if len(sys.argv) > 2 else 16)
    else:
        compare()
