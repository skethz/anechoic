"""PILOT sensitivity: can any plausible setting of ReAIM's unreported K2000 knobs reach p~0.47 / 0.80?"""
import json, time, itertools
from pathlib import Path
import numpy as np
import pilot as P, statica_repro as sr
Jf, edges = sr.load_graph(); J = Jf.astype(np.float32); sumw = float(edges[2].sum())
out = {}
for (fname, F), (t0, t1), kset, iters in itertools.product(
        [("max", np.max), ("min", np.min)], [(1.0, 0.1), (0.1, 0.01), (0.3, 0.003)],
        [(1, 2, 4, 8), (4, 8, 16, 32), (16, 32, 64, 128)], [1024, 4096]):
    rng = np.random.default_rng(abs(hash((fname, t0, kset, iters))) % 2**32)
    t = time.perf_counter()
    best, _ = P.asa(J, 128, iters, kset, 32, 96, rng, t0, t1, F)
    cut = (sumw - best) / 2; k = int((cut >= 33000).sum())
    key = f"F{fname}_T{t0}-{t1}_k{'-'.join(map(str,kset))}_it{iters}"
    out[key] = dict(F=fname, T=(t0, t1), kset=kset, iters=iters, successes=k, trials=128, p=k/128, mean_best_cut=float(cut.mean()), max_cut=float(cut.max()), sec=time.perf_counter()-t)
    print(f"{key:40s} p={k/128:.3f} mean={cut.mean():8.1f} max={cut.max():7.0f}", flush=True)
Path("sensitivity_results.json").write_text(json.dumps(out, indent=1) + "\n")
