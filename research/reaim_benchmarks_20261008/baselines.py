"""Descriptive baselines for reading normalized quality (not a method): expected normalized quality of a uniformly random
solution. MCP: random spins (cut ~ m/2); GPP: uniformly random exact bisection; TSP: uniformly random tour.
10,000 samples each, seed 20261008. Writes data/random_baselines.json."""
import json
import sys
from pathlib import Path
sys.dont_write_bytecode = True
HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import numpy as np  # noqa: E402
import problems as PB  # noqa: E402
import tsplib  # noqa: E402

rng = np.random.default_rng(20261008)
out = {}
bkv = json.loads((HERE.parent / 'gset_20261007' / 'bkv.json').read_text())['instances']
ref = json.loads((HERE / 'data' / 'gpp_reference.json').read_text())['instances']
for g in range(1, 21):
    P = PB.mcp(g, bkv[f'G{g}'])
    s = rng.choice(np.array([-1.0, 1.0], np.float32), size=(10000, P.N))
    out[f'mcp/G{g}'] = float(P.evaluate(s)['quality'].mean())
for g in (1, 2, 3, 4, 5, 14, 15, 16, 17):
    P = PB.gpp(g, 4, ref=ref[f'G{g}']['R'], target=ref[f'G{g}']['target'])
    s = np.ones((10000, P.N), np.float32)
    for r in range(10000):
        s[r, rng.permutation(P.N)[:P.N // 2]] = -1
    out[f'gpp/G{g}'] = float(P.evaluate(s)['quality'].mean())
for t in tsplib.INSTANCES:
    W = tsplib.load(t)['W']; n = len(W)
    L = np.array([W[p, np.roll(p, -1)].sum() for p in (rng.permutation(n) for _ in range(10000))])
    out[f'tsp/{t}'] = float((tsplib.OPT[t] / L).mean())
(HERE / 'data' / 'random_baselines.json').write_text(json.dumps(out, indent=1) + '\n')
print(json.dumps(out, indent=1))
