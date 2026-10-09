import sys, json, numpy as np
sys.path.insert(0, str(__import__('pathlib').Path(__file__).resolve().parent / '../theory_ideas_20261003'))
import sweeps_abc as m
import dmft
S = 200; Tt = m.sched(30.0, S) / np.sqrt(2000); qt = 4.0 / np.sqrt(2000)   # fixed rescaled schedule
pts = (25, 50, 100, 150, 200); out = {}
for N in (1000, 2000, 4000):
    vals = []
    for inst in range(4 if N < 4000 else 2):
        rng = np.random.default_rng(1000 * N + inst)
        W = np.triu(rng.choice([-1.0, 1.0], size=(N, N)).astype(np.float32), 1); J = W + W.T
        r = dmft.simulate(J, Tt * np.sqrt(N), np.full(S, qt * np.sqrt(N)), B=256 if N < 4000 else 128, seed=inst)
        vals.append([r['su'][t] for t in pts])
    v = np.array(vals); out[N] = v.tolist()
    print(f"N={N}: mean over instances " + " ".join(f"{x:.4f}" for x in v.mean(0)) + "   instance sd " + " ".join(f"{x:.4f}" for x in v.std(0)), flush=True)
json.dump(out, open('finite_n.json', 'w'))
