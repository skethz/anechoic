import sys, time, numpy as np
sys.path.insert(0, str(__import__('pathlib').Path(__file__).resolve().parent / '../theory_ideas_20261003'))
import sweeps_abc as m
import dmft
J, sumw = m.load()
S = 200; T = m.sched(30.0, S); q = np.full(S, 4.0)
r = dmft.simulate(J, T, q, B=1024, seed=2)
pts = (25, 50, 100, 150, 200)
print("sim B=1024:            " + " ".join(f"{r['su'][t]:.4f}" for t in pts))
for M in (50000,):
    for seed in (1, 2, 3):
        t0 = time.time(); d = dmft.solve_iter(T, q, M=M, iters=30, seed=seed)
        h = d['hist']
        print(f"iter M={M} s{seed} ({time.time()-t0:5.1f}s): " + " ".join(f"{d['su'][t]:.4f}" for t in pts)
              + f"  | last dC {h[-1]['dC']:.3f} dG {h[-1]['dG']:.3f}; su_final by iter: " + " ".join(f"{x['su_final']:.3f}" for x in h[::5]))
