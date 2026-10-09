"""The real Neal (dwave-samplers 1.2.0 SimulatedAnnealingSampler, the implementation behind `neal` 0.6.0) on a Max-Cut
instance, at Neal's defaults except num_sweeps and num_reads. Runs in the separate scratch venv venv_neal (numpy 1.26.4,
dimod 0.12.17, dwave-samplers 1.2.0); does not import any of this study's code.
Usage: python3 neal_ref.py <G-set file | WK2000_1.rud> <num_sweeps> <num_reads> <seed>  -> one JSON line on stdout."""
import json
import sys
import time

import numpy as np
import dimod
from dwave.samplers import SimulatedAnnealingSampler


def edges(path):
    lines = open(path).read().split('\n')
    n, m = (int(x) for x in lines[0].split())
    a = np.array([ln.split() for ln in lines[1:] if ln.strip()], dtype=np.int64)
    assert len(a) == m
    return n, a[:, 0] - 1, a[:, 1] - 1, a[:, 2]


def main():
    path, S, B, seed = sys.argv[1], int(sys.argv[2]), int(sys.argv[3]), int(sys.argv[4])
    n, i, j, w = edges(path)
    # maximize cut = sum w (1 - s_i s_j)/2  <=>  minimize sum w s_i s_j: Ising J_ij = w_ij, h = 0
    bqm = dimod.BinaryQuadraticModel.from_ising({k: 0.0 for k in range(n)}, {(int(a), int(b)): float(c) for a, b, c in zip(i, j, w)})
    t0 = time.time()
    ss = SimulatedAnnealingSampler().sample(bqm, num_reads=B, num_sweeps=S, seed=seed)
    order = list(ss.variables)
    X = np.asarray(ss.record.sample)[:, np.argsort(order)]          # columns in variable order 0..n-1
    cut = ((X[:, i] != X[:, j]) * w[None, :]).sum(1)
    print(json.dumps(dict(path=path, n=n, num_sweeps=S, num_reads=B, seed=seed, beta_range=list(map(float, ss.info['beta_range'])),
                          beta_schedule_type=ss.info.get('beta_schedule_type'), cuts=cut.tolist(), sec=time.time() - t0)))


if __name__ == '__main__':
    main()
