"""Verify the downloaded G-set files and write manifest.json (SHA-256, header checks, graph statistics).

Each file: first line "n m", then m lines "i j w" (1-based). Checks: edge count equals the header, weights in {-1, +1},
1 <= i, j <= n, no self-loops, no duplicate undirected edges. Statistics used by the protocol:
  W      = sum of edge weights (cut(s) = (W - H(s)) / 2 with H(s) = sum_{i<j} w_ij s_i s_j)
  dbar   = (1/N) sum_i sum_j J_ij^2 = 2m/N for unit weights (mean squared row norm = mean degree)
  sigma  = sqrt(dbar): RMS local field h_i = sum_j J_ij s_j under uniformly random spins
"""
import hashlib
import json
import math
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parent
DATA = ROOT / 'data'
URL = 'https://web.stanford.edu/~yyye/yyye/Gset/'
IDS = list(range(1, 48)) + [51, 52, 53, 54]
CLASSES = [  # (name, ids, N, structure, weights) as documented by the G-set generator (rudy) parameters
    ('R800+', range(1, 6), 800, 'random, density 6%', '+1'),
    ('R800+-', range(6, 11), 800, 'random, density 6% (graphs of G1-G5)', '+-1'),
    ('T800+-', range(11, 14), 800, 'toroidal 2D grid', '+-1'),
    ('P800+', range(14, 18), 800, 'planar-like (union of two planar graphs)', '+1'),
    ('P800+-', range(18, 22), 800, 'planar-like (graphs of G14-G17)', '+-1'),
    ('R2000+', range(22, 27), 2000, 'random, density 1%', '+1'),
    ('R2000+-', range(27, 32), 2000, 'random, density 1% (graphs of G22-G26)', '+-1'),
    ('T2000+-', range(32, 35), 2000, 'toroidal 2D grid', '+-1'),
    ('P2000+', range(35, 39), 2000, 'planar-like', '+1'),
    ('P2000+-', range(39, 43), 2000, 'planar-like (graphs of G35-G38)', '+-1'),
    ('R1000+', range(43, 48), 1000, 'random, density 2%', '+1'),
    ('P1000+', range(51, 55), 1000, 'planar-like', '+1'),
]


def cls_of(g):
    for name, ids, *_ in CLASSES:
        if g in ids:
            return name
    raise KeyError(g)


def parse(path):
    lines = path.read_text().split('\n')
    n, m = (int(x) for x in lines[0].split())
    rows = [ln.split() for ln in lines[1:] if ln.strip()]
    arr = np.array(rows, dtype=np.int64)
    return n, m, arr


def main():
    out = {'source_url': URL, 'download_utc': (ROOT / 'download_time_utc.txt').read_text().strip(),
           'note': 'G23 is listed on the server index with last-modified 2023-05-05 19:19; all other files 2003-09-11 12:16 '
                   '(see sources/gset_index.html).', 'instances': {}}
    ok_all = True
    for g in IDS:
        p = DATA / f'G{g}'
        raw = p.read_bytes()
        n, m, a = parse(p)
        i, j, w = a[:, 0], a[:, 1], a[:, 2]
        lo, hi = np.minimum(i, j), np.maximum(i, j)
        key = lo * (n + 1) + hi
        checks = dict(edge_count_matches_header=bool(len(a) == m), weights_pm1=bool(np.isin(w, (-1, 1)).all()),
                      indices_in_range=bool((i >= 1).all() and (j >= 1).all() and (i <= n).all() and (j <= n).all()),
                      no_self_loops=bool((i != j).all()), no_duplicate_edges=bool(len(np.unique(key)) == len(key)),
                      N_le_2048=bool(n <= 2048))
        deg = np.bincount(np.concatenate([i, j]) - 1, minlength=n)
        dbar = 2.0 * m / n
        rec = dict(file=f'data/G{g}', url=URL + f'G{g}', bytes=len(raw), sha256=hashlib.sha256(raw).hexdigest(),
                   N=n, m=m, n_pos=int((w > 0).sum()), n_neg=int((w < 0).sum()), W=int(w.sum()), cls=cls_of(g),
                   deg_min=int(deg.min()), deg_max=int(deg.max()), deg_mean=float(deg.mean()), isolated=int((deg == 0).sum()),
                   dbar=dbar, sigma=math.sqrt(dbar), checks=checks)
        ok = all(checks.values())
        ok_all &= ok
        out['instances'][f'G{g}'] = rec
        print(f"G{g:<3d} {rec['cls']:8s} N={n:5d} m={m:6d} +{rec['n_pos']:6d} -{rec['n_neg']:6d} W={rec['W']:6d} "
              f"deg {rec['deg_min']}-{rec['deg_max']} mean {dbar:6.2f} sigma {rec['sigma']:.3f} isolated {rec['isolated']} "
              f"{'OK' if ok else 'FAIL ' + str(checks)} {rec['sha256'][:16]}")
    out['all_checks_pass'] = bool(ok_all)
    out['classes'] = [dict(name=c[0], ids=[f'G{x}' for x in c[1]], N=c[2], structure=c[3], weights=c[4]) for c in CLASSES]
    (ROOT / 'manifest.json').write_text(json.dumps(out, indent=1) + '\n')
    print('ALL CHECKS PASS' if ok_all else 'SOME CHECKS FAILED')
    return 0 if ok_all else 1


if __name__ == '__main__':
    sys.exit(main())
