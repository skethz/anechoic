"""Ising instances for the ReAIM benchmark suite (MCP, GPP, TSP), in the engine convention

    E(s) = -sum_{i<j} J_ij s_i s_j - sum_i b_i s_i,      local field h_i = sum_{j != i} J_ij s_j + b_i,   dE_flip(i) = 2 s_i h_i.

J is stored as a sparse integer part (CSR, float32 values) plus an optional uniform part gamma on every off-diagonal pair:
J_ij = Js_ij + gamma (i != j). The uniform part carries the dense balance penalty of graph partitioning exactly, so every
field is an exact integer computed as h_i = (Js s + b)_i + gamma * (M - s_i), M = sum_j s_j. All values are integers below
2^24, so float32 arithmetic is exact and equal to a dense-matrix computation.

Formulations (Lucas 2014, Front. Phys. 2:5, Secs. 2.2 and 7.1), integer-scaled:
  MCP  J = -w (G-set), b = 0, as research/gset_20261007/engine.py. cut = (W - sum_{i<j} w_ij s_i s_j)/2, computed exactly.
  GPP  H = A (sum_i s_i)^2 + B sum_{(ij) in E} w_ij (1 - s_i s_j)/2, B = 1. Times 2: Js_ij = w_ij on edges, gamma = -P with
       P = 4A, b = 0 (N even). E(s) = 2 cut(s) + (P/2) M^2 + const. Feasible iff M = 0 (exact bisection).
  TSP  H = A sum_v (1 - sum_j x_vj)^2 + A sum_j (1 - sum_v x_vj)^2 + B sum_j sum_{u != v} W_uv x_uj x_v,j+1, B = 1,
       x_vj = (1 + s_vj)/2, spin index v*n + j (city v at tour position j, positions cyclic). Times 4:
       J = -P (P = 2A) on same-city and same-position pairs, J = -W_uv on (u, j), (v, j +- 1), u != v;
       b_vj = -2P(n - 2) - 2 D_v with D_v = sum_u W_uv. E(s) = 4 H(x) + const. Feasible iff x is a permutation matrix.
"""
import hashlib
import json
import math
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parent
GSET_DATA = ROOT.parent / 'gset_20261007' / 'data'
DN_K2000 = 1999.0 / 2000.0


class Problem:
    """Integer Ising instance with CSR sparse couplings, a uniform off-diagonal coupling gamma and a bias b."""

    def __init__(self, name, kind, N, rows, cols, vals, gamma, b, sigma_T, meta):
        self.name, self.kind, self.N = name, kind, int(N)
        rows = np.asarray(rows, np.int64); cols = np.asarray(cols, np.int64); vals = np.asarray(vals, np.float64)
        assert np.all(rows != cols)
        keep = vals != 0
        rows, cols, vals = rows[keep], cols[keep], vals[keep]
        order = np.lexsort((cols, rows))
        rows, cols, vals = rows[order], cols[order], vals[order]
        if len(rows) > 1:
            dup = (np.diff(rows) == 0) & (np.diff(cols) == 0)
            assert not dup.any(), 'duplicate coupling'
        self.indptr = np.zeros(self.N + 1, np.int64)
        np.add.at(self.indptr, rows + 1, 1)
        self.indptr = np.cumsum(self.indptr)
        self.indices = cols.astype(np.int64)
        self.data = vals.astype(np.float32)
        assert np.all(self.data.astype(np.float64) == vals), 'non-float32-exact coupling'
        self.gamma = float(gamma)
        self.b = np.asarray(b, np.float64).astype(np.float32)
        assert np.all(self.b.astype(np.float64) == np.asarray(b, np.float64))
        self.sigma_T = float(sigma_T)
        self.meta = meta
        self._stats = None

    # ---------------------------------------------------------------- matrix access (small instances / export only)
    def dense(self):
        J = np.full((self.N, self.N), self.gamma, np.float64)
        np.fill_diagonal(J, 0.0)
        r = np.repeat(np.arange(self.N), np.diff(self.indptr))
        J[r, self.indices] += self.data
        return J

    def check_symmetric(self):
        r = np.repeat(np.arange(self.N), np.diff(self.indptr))
        a = dict(zip(zip(r.tolist(), self.indices.tolist()), self.data.tolist()))
        return all(a.get((j, i)) == v for (i, j), v in a.items())

    # ---------------------------------------------------------------- energies and fields (exact, float64)
    def field(self, s):
        s = np.asarray(s, np.float64)
        B = s.shape[0]
        r = np.repeat(np.arange(self.N), np.diff(self.indptr))
        h = np.zeros((B, self.N))
        np.add.at(h.T, r, (self.data.astype(np.float64)[:, None] * s[:, self.indices].T))
        h += self.gamma * (s.sum(1, keepdims=True) - s) + self.b.astype(np.float64)[None, :]
        return h

    def energy(self, s):
        s = np.asarray(s, np.float64)
        h = self.field(s)
        return -0.5 * np.einsum('bi,bi->b', s, h - self.b[None, :]) - s @ self.b.astype(np.float64)

    # ---------------------------------------------------------------- statistics used for scaling and reporting
    def stats(self):
        if self._stats is None:
            N = self.N; g = self.gamma
            r = np.repeat(np.arange(N), np.diff(self.indptr))
            v = self.data.astype(np.float64)
            # row sums of J_ij^2 with J_ij = Js_ij + gamma off-diagonal
            rowsq = np.full(N, g * g * (N - 1))
            np.add.at(rowsq, r, (v + g) ** 2 - g * g)
            rowabs = np.full(N, abs(g) * (N - 1))
            np.add.at(rowabs, r, np.abs(v + g) - abs(g))
            levels = set(np.unique(np.abs(v + g)).tolist())
            if g != 0 and (N * (N - 1) - len(v)) > 0:
                levels.add(abs(g))
            levels.discard(0.0)
            maxabs = max(levels) if levels else 0.0
            nnz_offdiag = len(v) if g == 0 else N * (N - 1) - int(np.sum(v + g == 0))
            self._stats = dict(N=N, nnz_offdiag=int(nnz_offdiag), density=nnz_offdiag / (N * (N - 1)),
                               maxabs_J=float(maxabs), distinct_absJ_levels=len(levels),
                               absJ_levels=sorted(levels)[:64], mean_rowsq=float(rowsq.mean()),
                               max_rowabs=float(rowabs.max()), max_abs_b=float(np.abs(self.b).max()),
                               distinct_b=int(len(np.unique(self.b))),
                               field_bound=float((rowabs + np.abs(self.b.astype(np.float64))).max()),
                               sigma_T=self.sigma_T)
            st = self._stats
            st['bits_J_signed'] = int(1 + math.ceil(math.log2(maxabs + 1))) if maxabs > 0 else 1
            st['bits_b_signed'] = int(1 + math.ceil(math.log2(st['max_abs_b'] + 1))) if st['max_abs_b'] > 0 else 1
            st['bits_field_signed'] = int(1 + math.ceil(math.log2(st['field_bound'] + 1)))
            st['sigmaJ_goto'] = math.sqrt(rowsq.sum() / (N * (N - 1)))
            st['lam_factor'] = (st['mean_rowsq'] / N) / DN_K2000
            st['kappa_factor'] = st['mean_rowsq'] / self.sigma_T ** 2
        return self._stats

    # ---------------------------------------------------------------- objective (true, unquantized problem)
    def evaluate(self, s):
        """Per run: dict of arrays. 'feasible' (bool), 'value' (cut / balanced cut / tour length, -1 if infeasible),
        'quality' (normalized objective in [0, ~1], 0 if infeasible), 'success' (target reached), 'opt' (optimum or
        best-known reached)."""
        s = np.asarray(s)
        m = self.meta
        if self.kind == 'mcp':
            si = s[:, m['ei']]; sj = s[:, m['ej']]
            cut = ((si != sj).astype(np.int64) * m['w'][None, :]).sum(1)
            feas = np.ones(len(s), bool)
            return dict(feasible=feas, value=cut, quality=cut / m['ref'], success=cut >= m['target'],
                        opt=cut >= m['ref'])
        if self.kind == 'gpp':
            si = s[:, m['ei']]; sj = s[:, m['ej']]
            cut = ((si != sj).astype(np.int64) * m['w'][None, :]).sum(1)
            M = np.rint(s.astype(np.float64).sum(1)).astype(np.int64)
            feas = M == 0
            if m.get('ref') is None:      # reference not yet known (reference computation, calibration)
                nan = np.full(len(s), np.nan)
                return dict(feasible=feas, value=np.where(feas, cut, -1), cut_any=cut, imbalance=M, quality=nan,
                            success=np.zeros(len(s), bool), opt=np.zeros(len(s), bool))
            q = np.where(feas, m['ref'] / np.maximum(cut, 1), 0.0)
            return dict(feasible=feas, value=np.where(feas, cut, -1), cut_any=cut, imbalance=M, quality=q,
                        success=feas & (cut <= m['target']), opt=feas & (cut <= m['ref']))
        if self.kind == 'tsp':
            n = m['n']; W = m['W']
            x = (s.reshape(len(s), n, n) > 0)
            feas = np.all(x.sum(1) == 1, axis=1) & np.all(x.sum(2) == 1, axis=1)
            L = np.full(len(s), -1, np.int64)
            for r in np.nonzero(feas)[0]:
                tour = np.argmax(x[r], axis=0)             # city at each position
                L[r] = int(W[tour, np.roll(tour, -1)].sum())
            q = np.where(feas, m['ref'] / np.maximum(L, 1), 0.0)
            return dict(feasible=feas, value=L, quality=q, success=feas & (L <= m['target']), opt=feas & (L <= m['ref']),
                        rows_ok=np.all(x.sum(2) == 1, axis=1), cols_ok=np.all(x.sum(1) == 1, axis=1),
                        n_on=x.reshape(len(s), -1).sum(1))
        raise ValueError(self.kind)

    # ---------------------------------------------------------------- K-bit quantization (precision study, export)
    def quantize(self, K):
        """Couplings to signed K-bit integers |J| <= 2^(K-1) - 1: alpha = (2^(K-1) - 1) / max|J|, J' = round(alpha J),
        b' = round(alpha b), rounding half away from zero, applied entrywise to the full J (sparse part + gamma).
        K = None returns self (full precision, alpha = 1). The objective stays on the true problem (meta unchanged)."""
        if K is None:
            return self, 1.0
        st = self.stats()
        alpha = (2 ** (K - 1) - 1) / st['maxabs_J']
        rnd = lambda x: np.sign(x) * np.floor(np.abs(x) + 0.5)  # noqa: E731
        g2 = float(rnd(alpha * self.gamma)) if self.gamma != 0 else 0.0
        r = np.repeat(np.arange(self.N), np.diff(self.indptr))
        full = self.data.astype(np.float64) + self.gamma
        v2 = rnd(alpha * full) - g2
        b2 = rnd(alpha * self.b.astype(np.float64))
        q = Problem(f'{self.name}_K{K}', self.kind, self.N, r, self.indices, v2, g2, b2, self.sigma_T * alpha, self.meta)
        q.K = K
        q.alpha = alpha
        assert q.stats()['maxabs_J'] <= 2 ** (K - 1) - 1
        return q, alpha

    def jint8(self):
        """Dense int8 matrix (zero diagonal) for the SCAJINT8 export; requires |J| <= 127."""
        J = self.dense()
        assert np.abs(J).max() <= 127 and np.all(J == np.round(J))
        return J.astype(np.int8)


def write_jint8(path, J8):
    n = J8.shape[0]
    with open(path, 'wb') as f:
        f.write(b'SCAJINT8'); f.write(np.uint32(n).tobytes()); f.write(np.uint32(0).tobytes())
        f.write(np.ascontiguousarray(J8, dtype=np.int8).tobytes())


def write_bias(path, b):
    np.ascontiguousarray(np.asarray(b), dtype='<i4').tofile(path)


def sha256(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


# ==================================================================== builders
def _gset_edges(g):
    lines = (GSET_DATA / f'G{g}').read_text().split('\n')
    n, m = (int(x) for x in lines[0].split())
    a = np.array([ln.split() for ln in lines[1:] if ln.strip()], dtype=np.int64)
    assert len(a) == m
    return n, a[:, 0] - 1, a[:, 1] - 1, a[:, 2]


def mcp(g, bkv):
    n, i, j, w = _gset_edges(g)
    rows = np.concatenate([i, j]); cols = np.concatenate([j, i]); vals = -np.concatenate([w, w]).astype(np.float64)
    dbar = 2.0 * len(i) / n
    meta = dict(ei=i, ej=j, w=w, W=int(w.sum()), ref=int(bkv['BKV']), target=int(bkv['target']), m=len(i))
    return Problem(f'G{g}', 'mcp', n, rows, cols, vals, 0.0, np.zeros(n), math.sqrt(dbar), meta)


def gpp_from_edges(name, n, i, j, w, P, ref=None, target=None):
    """Graph bisection, Lucas Sec. 2.2 with B = 1, A = P/4 (integer form x2): J = w[E] - P, b = 0."""
    assert n % 2 == 0
    rows = np.concatenate([i, j]); cols = np.concatenate([j, i]); vals = np.concatenate([w, w]).astype(np.float64)
    dbar = float((w.astype(np.float64) ** 2).sum() * 2 / n)
    meta = dict(ei=i, ej=j, w=w, m=len(i), P=P, A=P / 4.0, ref=ref, target=target, dbar=dbar)
    return Problem(name, 'gpp', n, rows, cols, vals, -float(P), np.zeros(n), math.sqrt(dbar), meta)


def gpp(g, P, ref=None, target=None):
    n, i, j, w = _gset_edges(g)
    return gpp_from_edges(f'G{g}', n, i, j, w, P, ref, target)


def tsp_from_matrix(name, W, A, ref=None, target=None):
    """TSP, Lucas Sec. 7.1 with B = 1 (integer form x4): P = 2A. Spin index v*n + j."""
    W = np.asarray(W, np.int64); n = W.shape[0]
    assert float(2 * A) == round(2 * A), 'A must be a multiple of 1/2'
    P = int(round(2 * A))
    idx = lambda v, j: v * n + (j % n)  # noqa: E731
    rows, cols, vals = [], [], []
    for v in range(n):
        for j in range(n):
            a = idx(v, j)
            for k in range(n):               # same city, other positions
                if k != j:
                    rows.append(a); cols.append(idx(v, k)); vals.append(-P)
            for u in range(n):               # same position, other cities
                if u != v:
                    rows.append(a); cols.append(idx(u, j)); vals.append(-P)
            for u in range(n):               # adjacent positions, other cities
                if u != v:
                    for dj in (1, -1):
                        rows.append(a); cols.append(idx(u, j + dj)); vals.append(-float(W[v, u]))
    D = W.sum(1)
    b = np.array([-2.0 * P * (n - 2) - 2.0 * D[v] for v in range(n) for j in range(n)])
    meta = dict(n=n, W=W, A=A, P=P, ref=ref, target=target)
    return Problem(name, 'tsp', n * n, rows, cols, vals, 0.0, b, float(A), meta)


def tsp(name, A_rule='MQC'):
    import tsplib
    d = tsplib.load(name)
    W = d['W']
    A = float(W.max()) if A_rule == 'MQC' else float(A_rule)
    opt = tsplib.OPT[name]
    return tsp_from_matrix(name, W, A, ref=opt, target=int(math.floor(1.01 * opt)))


def tsp_constant(W, A):
    """Constant c with E(s) = 4 H(x) + c (used only by the brute-force verification)."""
    n = W.shape[0]
    # group penalty constant per group: 1 - n/2 + C(n,2)/2 ; 2n groups ; times A ; distance constant sum_j sum_{u!=v} W/4
    const_H = A * 2 * n * (1 - n / 2 + n * (n - 1) / 4) + n * (W.sum()) / 4.0
    return -4.0 * const_H
