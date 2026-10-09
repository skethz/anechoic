"""PILOT: larger flip caps with the paper's Max-Cut settings (F=max, T 1->0.1)."""
import json, time, itertools
from pathlib import Path
import numpy as np
import pilot as P, statica_repro as sr
Jf, edges = sr.load_graph(); J = Jf.astype(np.float32); sumw = float(edges[2].sum())
out = {}
for kset, iters in itertools.product([(64, 128, 256, 512), (128, 256, 512, 1024), (256, 512, 1024, 2000)], [512, 1024, 2048, 4096]):
    rng = np.random.default_rng(abs(hash((kset, iters, 'big'))) % 2**32)
    t = time.perf_counter()
    best, _ = P.asa(J, 128, iters, kset, 32, 96, rng, 1.0, 0.1, np.max)
    cut = (sumw - best) / 2; k = int((cut >= 33000).sum())
    key = f"Fmax_T1.0-0.1_k{'-'.join(map(str,kset))}_it{iters}"
    out[key] = dict(kset=kset, iters=iters, successes=k, trials=128, p=k/128, p_ci95=sr.clopper_pearson(k,128), mean_best_cut=float(cut.mean()), max_cut=float(cut.max()), sec=time.perf_counter()-t)
    print(f"{key:44s} p={k/128:.3f} mean={cut.mean():8.1f} max={cut.max():7.0f} ({time.perf_counter()-t:.0f}s)", flush=True)
Path("sensitivity_bigk_results.json").write_text(json.dumps(out, indent=1) + "\n")
