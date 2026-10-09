"""Ancillary-spin reduction of an Ising problem with local fields (Goto et al. 2021, main text and section S1: "introducing
an ancillary spin reduces the Ising problem to the one without local fields"): spin 0' = index N couples to every spin i
with J'_{iN} = b_i; the new problem has no fields. A solution s' maps back by s_i = s'_i s'_N (global-flip symmetry).
Used for the SB methods (aSB, bSB, dSB) on TSP, the only problem here with fields."""
import numpy as np

import problems as PB


def with_ancilla(P):
    N = P.N
    r = np.repeat(np.arange(N), np.diff(P.indptr))
    nz = np.nonzero(P.b)[0]
    rows = np.concatenate([r, nz, np.full(len(nz), N)])
    cols = np.concatenate([P.indices, np.full(len(nz), N), nz])
    vals = np.concatenate([P.data.astype(np.float64), P.b[nz].astype(np.float64), P.b[nz].astype(np.float64)])
    assert P.gamma == 0.0, 'ancilla construction assumes no uniform coupling (TSP)'
    Q = PB.Problem(P.name + '+anc', P.kind, N + 1, rows, cols, vals, 0.0, np.zeros(N + 1), P.sigma_T, P.meta)
    return Q


def decode(s_anc):
    s = np.asarray(s_anc)
    return (s[:, :-1] * s[:, -1:]).astype(np.float32)
