import sys, json, numpy as np
sys.path.insert(0, str(__import__('pathlib').Path(__file__).resolve().parent / '../theory_ideas_20261003'))
import sweeps_abc as m
import dmft
S = 560; Tt = m.sched(30.0, S) / np.sqrt(2000); qt = 4.0 / np.sqrt(2000); pts = (140, 280, 420, 560); out = {}
for N in (1000, 2000, 4000, 8000):
    vals = []; fl = []
    for inst in range(3 if N <= 2000 else (2 if N == 4000 else 1)):
        rng = np.random.default_rng(7000 * N + inst)
        W = np.triu(rng.choice([-1.0, 1.0], size=(N, N)).astype(np.float32), 1); J = W + W.T
        r = dmft.simulate(J, Tt * np.sqrt(N), np.full(S, qt * np.sqrt(N)), lam=np.full(S, 0.7), B=128 if N <= 4000 else 48, seed=inst)
        vals.append([r['su'][t] for t in pts]); fl.append(r['flips'].sum())
    v = np.array(vals); out[N] = dict(su=v.tolist(), flips_per_spin=float(np.mean(fl)))
    print(f"N={N}: su[140,280,420,560] " + " ".join(f"{x:.4f}" for x in v.mean(0)) + f"  flips/spin {np.mean(fl):.1f}", flush=True)
json.dump(out, open('finite_n_lam.json', 'w'))
