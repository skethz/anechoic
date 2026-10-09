"""Synthetic graphs of the G-set classes (NOT G-set instances), used only for engine verification and the exploratory
pre-protocol calibration. Generators mimic the G-set classes: random G(n, m), 2D torus, planar-like (union of two random
planar triangulations), with +1 or +-1 weights."""
import math

import numpy as np

from engine import Graph


def _weights(rng, m, signed):
    return rng.choice(np.array([-1, 1]), size=m) if signed else np.ones(m, np.int64)


def random_graph(n, m, signed, seed):
    rng = np.random.default_rng(seed)
    iu = np.triu_indices(n, 1)
    pick = rng.choice(len(iu[0]), size=m, replace=False)
    i, j = iu[0][pick], iu[1][pick]
    return Graph(f'syn_rand_n{n}_m{m}_{"pm" if signed else "p"}_{seed}', n, (i, j, _weights(rng, m, signed)))


def torus(L1, L2, seed):
    rng = np.random.default_rng(seed)
    idx = np.arange(L1 * L2).reshape(L1, L2)
    i = np.concatenate([idx.ravel(), idx.ravel()])
    j = np.concatenate([np.roll(idx, -1, 0).ravel(), np.roll(idx, -1, 1).ravel()])
    return Graph(f'syn_torus_{L1}x{L2}_{seed}', L1 * L2, (i, j, _weights(rng, len(i), True)))


def planar_like(n, signed, seed):
    """Union of two planar graphs on the same vertex set: each is an L x L grid with one random diagonal per cell
    (a planar triangulation), with an independent random vertex labelling (scipy.spatial is unavailable here)."""
    rng = np.random.default_rng(seed)
    L = int(math.ceil(math.sqrt(n)))
    E = set()
    for _ in range(2):
        lab = rng.permutation(L * L)
        idx = lab.reshape(L, L)
        cand = []
        for r in range(L):
            for c in range(L):
                if c + 1 < L: cand.append((idx[r, c], idx[r, c + 1]))
                if r + 1 < L: cand.append((idx[r, c], idx[r + 1, c]))
                if r + 1 < L and c + 1 < L:
                    cand.append((idx[r, c], idx[r + 1, c + 1]) if rng.random() < 0.5 else (idx[r, c + 1], idx[r + 1, c]))
        for x, y in cand:
            if x < n and y < n and x != y:
                E.add((int(min(x, y)), int(max(x, y))))
    e = np.array(sorted(E))
    return Graph(f'syn_planar_n{n}_{"pm" if signed else "p"}_{seed}', n, (e[:, 0], e[:, 1], _weights(rng, len(e), signed)))
