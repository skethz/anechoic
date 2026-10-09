"""Reproducibility and speed check: re-run original pilot points with their ORIGINAL seeds and compare with results/S250.json
(no new data is selected from this; it only checks that this machine reproduces the laptop's numbers)."""
import json, sys, time
from pathlib import Path
import numpy as np
sys.path.insert(0, str(Path(__file__).resolve().parent))
from common import M, ROOT
import importlib.util
spec = importlib.util.spec_from_file_location('sweep_orig', ROOT / 'research/algorithm_compare_20261007/sweep.py')
sw = importlib.util.module_from_spec(spec); spec.loader.exec_module(sw)
orig = json.loads(Path(sys.argv[1]).read_text())
np.seterr(all='ignore')
J, sumw = M.m.load()
S = 250; si = 0
for mi, (name, grid) in enumerate(sw.grids(S)):
    gi = orig[name]['selected']
    cfg = grid[gi]
    t = time.time()
    E, cuts = M.run_method(J, sumw, cfg, 64, 100000 * mi + 80000 + 1000 * si + gi)
    dt = time.time() - t
    o = orig[name]['pilot'][gi]['mean_cut']
    print(f"{name:20s} gi={gi:2d} mean_cut={cuts.mean():.4f} orig={o:.4f} match={abs(cuts.mean()-o)<1e-9} sec={dt:.1f}", flush=True)
