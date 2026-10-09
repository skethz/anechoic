import sys, time, json, numpy as np
sys.path.insert(0, str(__import__('pathlib').Path(__file__).resolve().parent / '../theory_ideas_20261003'))
import sweeps_abc as m
import dmft
S = 560; T = m.sched(30.0, S); q = np.full(S, 4.0)
sim = {0.0: 1.4918, 0.7: 1.5004}
for lam in (0.0, 0.7):
    for seed in (1, 2):
        t0 = time.time()
        d = dmft.solve_iter(T, q, lam=np.full(S, lam), M=30000, iters=200, alpha=0.05, avg_last=60, seed=seed)
        tr = [x['su_final'] for x in d['hist']]
        print(f"lambda {lam} s{seed} ({time.time()-t0:.0f}s): final su {d['su'][-1]:.4f} (sim {sim[lam]}) su[280] {d['su'][280]:.4f} | "
              "iter-means 0-200 step 25: " + " ".join(f"{np.mean(tr[i:i+25]):.3f}" for i in range(0, 200, 25)), flush=True)
