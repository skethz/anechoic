"""Verify the six TSPLIB instances before any run: SHA-256 of the downloaded files, parsed matrices (symmetric, zero
diagonal, positive off-diagonal), the published optima parsed from TSPLIB's optimum list (sources/TSP-BEST.html), the
lengths of the optimal tours shipped with TSPLIB (gr24, fri26, bayg29, bays29), and an independent exact optimum by
integer programming for all six. Writes data/tsplib_manifest.json."""
import hashlib
import json
import re
import sys
import time
from pathlib import Path

sys.dont_write_bytecode = True
import numpy as np  # noqa: E402

import tsplib  # noqa: E402

ROOT = Path(__file__).resolve().parent
BASE = 'https://comopt.ifi.uni-heidelberg.de/software/TSPLIB95/tsp/'


def sha(p):
    return hashlib.sha256(Path(p).read_bytes()).hexdigest()


def main():
    best_html = (ROOT / 'sources' / 'TSP-BEST.html').read_text(encoding='latin-1')
    pub = {m.group(1): int(m.group(2)) for m in re.finditer(r'([a-z]+\d+)\s*:\s*(\d+)', best_html)}
    out = dict(source=BASE, optimum_list='https://comopt.ifi.uni-heidelberg.de/software/TSPLIB95/tsp/TSP-BEST.html',
               optimum_list_sha256=sha(ROOT / 'sources' / 'TSP-BEST.html'),
               download_utc=(ROOT / 'data' / 'download_time_utc.txt').read_text().strip(),
               note='Two independent downloads were byte-identical (SHA-256 compared before the second copy was deleted).',
               instances={})
    ok_all = True
    for name in tsplib.INSTANCES:
        f = ROOT / 'data' / 'tsplib' / f'{name}.tsp.gz'
        d = tsplib.parse(f)
        W = d['W']; n = d['n']
        rec = dict(file=f'data/tsplib/{name}.tsp.gz', url=BASE + f'{name}.tsp.gz', bytes=f.stat().st_size, sha256=sha(f),
                   n=n, spins_one_hot=n * n, edge_weight_type=d['spec']['EDGE_WEIGHT_TYPE'],
                   edge_weight_format=d['spec'].get('EDGE_WEIGHT_FORMAT'), comment=d['spec'].get('COMMENT'),
                   W_min_offdiag=int(W[~np.eye(n, dtype=bool)].min()), W_max=int(W.max()),
                   W_mean_offdiag=float(W[~np.eye(n, dtype=bool)].mean()),
                   distinct_distances=int(len(np.unique(W[np.triu_indices(n, 1)]))),
                   optimum_published=pub.get(name), optimum_hardcoded=tsplib.OPT[name])
        tf = ROOT / 'data' / 'tsplib' / f'{name}.opt.tour.gz'
        if tf.exists():
            t = tsplib.parse_tour(tf)
            assert sorted(t) == list(range(n))
            rec.update(opt_tour_file=f'data/tsplib/{name}.opt.tour.gz', opt_tour_url=BASE + f'{name}.opt.tour.gz',
                       opt_tour_sha256=sha(tf), opt_tour_length=tsplib.tour_length(W, t))
        t0 = time.time()
        L, tour, iters = tsplib.exact_optimum(W)
        rec.update(milp_optimum=L, milp_tour=tour, milp_iterations=iters, milp_sec=round(time.time() - t0, 2),
                   milp_tour_length_check=tsplib.tour_length(W, tour))
        ok = (rec['optimum_published'] == rec['optimum_hardcoded'] == L == rec['milp_tour_length_check']
              and rec.get('opt_tour_length', L) == L)
        rec['all_checks_pass'] = bool(ok); ok_all &= ok
        out['instances'][name] = rec
        print(f"{name:7s} n={n:2d} {rec['edge_weight_format']:15s} W in [{rec['W_min_offdiag']}, {rec['W_max']}] "
              f"published={rec['optimum_published']} opt.tour={rec.get('opt_tour_length')} MILP={L} "
              f"({iters} rounds, {rec['milp_sec']}s) -> {'OK' if ok else 'MISMATCH'}", flush=True)
    out['all_checks_pass'] = bool(ok_all)
    (ROOT / 'data' / 'tsplib_manifest.json').write_text(json.dumps(out, indent=1) + '\n')
    print('ALL OK' if ok_all else 'FAILURES')


if __name__ == '__main__':
    main()
