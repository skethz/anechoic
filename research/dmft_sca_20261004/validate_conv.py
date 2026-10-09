import sys, time, numpy as np
sys.path.insert(0, str(__import__('pathlib').Path(__file__).resolve().parent / '../theory_ideas_20261003'))
import sweeps_abc as m
import dmft
J, sumw = m.load()
S = 200; T = m.sched(30.0, S); q = np.full(S, 4.0)
r = dmft.simulate(J, T, q, B=1024, seed=2)
print("sim B=1024:      " + " ".join(f"{r['su'][t]:.4f}" for t in (25, 50, 100, 150, 200)))
for est in ('lr', 'rb'):
    for M in (40000, 160000):
        for seed in (1, 2):
            t0 = time.time(); d = dmft.solve(T, q, M=M, seed=seed, estimator=est)
            print(f"{est} M={M:6d} s{seed} ({time.time()-t0:4.1f}s): " + " ".join(f"{d['su'][t]:.4f}" for t in (25, 50, 100, 150, 200)))
