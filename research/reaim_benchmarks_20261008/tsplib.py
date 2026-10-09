"""TSPLIB parser (explicit and geographic formats, per Reinelt's TSPLIB 95 documentation, sources/tsp95.pdf) and an exact
optimum check by integer programming (scipy HiGHS, undirected edge variables, degree constraints, iterative subtour
elimination). Used to verify the published optima of gr17, gr21, gr24, fri26, bayg29 and bays29 before any run.

Distances follow TSPLIB 95 Sec. 2: EXPLICIT matrices are integers as given; GEO uses nint() of the coordinates' degrees and
minutes, RRR = 6378.388, and dij = (int)(RRR * acos(0.5*((1+q1)*q2 - (1-q1)*q3)) + 1.0); EUC_2D uses nint(sqrt(.)).
"""
import gzip
import math
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parent
DATA = ROOT / 'data' / 'tsplib'
INSTANCES = ('gr17', 'gr21', 'gr24', 'fri26', 'bayg29', 'bays29')
# published optima, TSPLIB "Optimal solutions for symmetric TSPs" (sources/TSP-BEST.html, parsed by verify_tsplib.py)
OPT = {'gr17': 2085, 'gr21': 2707, 'gr24': 1272, 'fri26': 937, 'bayg29': 1610, 'bays29': 2020}

KEYWORD_SECTIONS = ('NODE_COORD_SECTION', 'EDGE_WEIGHT_SECTION', 'DISPLAY_DATA_SECTION', 'TOUR_SECTION',
                    'FIXED_EDGES_SECTION', 'DEPOT_SECTION', 'DEMAND_SECTION', 'EOF')


def _read(path):
    p = Path(path)
    raw = gzip.decompress(p.read_bytes()) if p.suffix == '.gz' else p.read_bytes()
    return raw.decode('ascii')


def parse(path):
    """Returns dict(name, n, W (n x n int64, symmetric, zero diagonal), spec (header dict))."""
    lines = _read(path).splitlines()
    spec, k = {}, 0
    sections = {}
    while k < len(lines):
        ln = lines[k].strip(); k += 1
        if not ln:
            continue
        key = ln.split(':')[0].strip() if ':' in ln else ln.split()[0]
        if key in KEYWORD_SECTIONS:
            if key == 'EOF':
                break
            body = []
            while k < len(lines):
                t = lines[k].strip()
                tk = t.split(':')[0].strip() if ':' in t else (t.split()[0] if t else '')
                if tk in KEYWORD_SECTIONS:
                    break
                body.append(t); k += 1
            sections[key] = ' '.join(body).split()
        else:
            spec[key] = ln.split(':', 1)[1].strip()
    n = int(spec['DIMENSION'])
    ewt = spec['EDGE_WEIGHT_TYPE']
    W = np.zeros((n, n), np.int64)
    if ewt == 'EXPLICIT':
        fmt = spec['EDGE_WEIGHT_FORMAT']
        v = [int(x) for x in sections['EDGE_WEIGHT_SECTION']]
        it = iter(v)
        if fmt == 'FULL_MATRIX':
            assert len(v) == n * n, (len(v), n)
            W = np.array(v, np.int64).reshape(n, n)
        elif fmt == 'UPPER_ROW':            # row i: j = i+1..n-1
            assert len(v) == n * (n - 1) // 2
            for i in range(n):
                for j in range(i + 1, n):
                    W[i, j] = W[j, i] = next(it)
        elif fmt == 'LOWER_ROW':            # row i: j = 0..i-1
            assert len(v) == n * (n - 1) // 2
            for i in range(n):
                for j in range(i):
                    W[i, j] = W[j, i] = next(it)
        elif fmt == 'UPPER_DIAG_ROW':       # row i: j = i..n-1
            assert len(v) == n * (n + 1) // 2
            for i in range(n):
                for j in range(i, n):
                    W[i, j] = W[j, i] = next(it)
        elif fmt == 'LOWER_DIAG_ROW':       # row i: j = 0..i
            assert len(v) == n * (n + 1) // 2, (len(v), n)
            for i in range(n):
                for j in range(i + 1):
                    W[i, j] = W[j, i] = next(it)
        else:
            raise ValueError(fmt)
    elif ewt in ('GEO', 'EUC_2D', 'ATT', 'CEIL_2D'):
        c = np.array(sections['NODE_COORD_SECTION'], float).reshape(n, 3)[:, 1:]
        for i in range(n):
            for j in range(n):
                if i != j:
                    W[i, j] = _dist(ewt, c[i], c[j])
    else:
        raise ValueError(ewt)
    assert np.all(W == W.T), 'asymmetric'
    assert np.all(np.diag(W) == 0), 'nonzero diagonal'
    assert np.all(W[~np.eye(n, dtype=bool)] > 0), 'nonpositive off-diagonal distance'
    return dict(name=spec['NAME'], n=n, W=W, spec=spec)


def _nint(x):
    return int(x + 0.5)


def _geo_rad(x):
    # TSPLIB FAQ (sources/TSPFAQ.html): PI = 3.141592; deg = (int) x; min = x - deg; rad = PI*(deg + 5.0*min/3.0)/180.0
    deg = int(x)
    mn = x - deg
    return 3.141592 * (deg + 5.0 * mn / 3.0) / 180.0


def _dist(ewt, a, b):
    if ewt == 'EUC_2D':
        return _nint(math.hypot(a[0] - b[0], a[1] - b[1]))
    if ewt == 'CEIL_2D':
        return int(math.ceil(math.hypot(a[0] - b[0], a[1] - b[1])))
    if ewt == 'ATT':
        r = math.sqrt(((a[0] - b[0]) ** 2 + (a[1] - b[1]) ** 2) / 10.0); t = _nint(r)
        return t + 1 if t < r else t
    if ewt == 'GEO':
        RRR = 6378.388
        lat_i, lon_i, lat_j, lon_j = _geo_rad(a[0]), _geo_rad(a[1]), _geo_rad(b[0]), _geo_rad(b[1])
        q1 = math.cos(lon_i - lon_j); q2 = math.cos(lat_i - lat_j); q3 = math.cos(lat_i + lat_j)
        return int(RRR * math.acos(0.5 * ((1.0 + q1) * q2 - (1.0 - q1) * q3)) + 1.0)
    raise ValueError(ewt)


def parse_tour(path):
    lines = _read(path).split()
    k = lines.index('TOUR_SECTION')
    t = []
    for x in lines[k + 1:]:
        v = int(x)
        if v == -1:
            break
        t.append(v - 1)
    return t


def tour_length(W, tour):
    t = list(tour)
    return int(sum(W[t[i], t[(i + 1) % len(t)]] for i in range(len(t))))


def exact_optimum(W, time_limit=600.0):
    """Exact symmetric TSP optimum by MILP with lazy subtour elimination (scipy.optimize.milp / HiGHS).
    Returns (length, tour). Iterates: solve with degree constraints + collected subtour cuts; if the solution has
    several cycles, add sum_{e in delta(S)} x_e >= 2 for every cycle S; repeat until a single Hamiltonian cycle."""
    from scipy.optimize import milp, LinearConstraint, Bounds
    from scipy.sparse import lil_matrix
    n = W.shape[0]
    E = [(i, j) for i in range(n) for j in range(i + 1, n)]
    m = len(E); c = np.array([W[i, j] for i, j in E], float)
    A = lil_matrix((n, m))
    for e, (i, j) in enumerate(E):
        A[i, e] = 1; A[j, e] = 1
    cons = [LinearConstraint(A.tocsr(), 2, 2)]
    cuts = []
    for it in range(500):
        allc = cons + ([LinearConstraint(np.array(cuts), 2, np.inf)] if cuts else [])
        res = milp(c, constraints=allc, integrality=np.ones(m), bounds=Bounds(0, 1),
                   options=dict(time_limit=time_limit, mip_rel_gap=0.0))
        assert res.status == 0, res.message
        x = np.round(res.x).astype(int)
        adj = [[] for _ in range(n)]
        for e, (i, j) in enumerate(E):
            if x[e]:
                adj[i].append(j); adj[j].append(i)
        seen = [False] * n; comps = []
        for s in range(n):
            if seen[s]:
                continue
            comp = []; st = [s]; seen[s] = True
            while st:
                u = st.pop(); comp.append(u)
                for v in adj[u]:
                    if not seen[v]:
                        seen[v] = True; st.append(v)
            comps.append(comp)
        if len(comps) == 1:
            tour = [0]; prev = -1
            while len(tour) < n:
                u = tour[-1]; nxt = [v for v in adj[u] if v != prev][0]
                prev = u; tour.append(nxt)
            return int(round(res.fun)), tour, it + 1
        for comp in comps:
            S = set(comp); row = np.zeros(m)
            for e, (i, j) in enumerate(E):
                if (i in S) != (j in S):
                    row[e] = 1
            cuts.append(row)
    raise RuntimeError('subtour loop did not converge')


def load(name):
    return parse(DATA / f'{name}.tsp.gz')
