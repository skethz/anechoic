import sys, time, numpy as np
sys.path.insert(0, str(__import__('pathlib').Path(__file__).resolve().parent / '../theory_ideas_20261003'))
import sweeps_abc as m
import dmft
J, sumw = m.load()
S = 200; T = m.sched(30.0, S); q = np.full(S, 4.0)
r = dmft.simulate(J, T, q, B=512, seed=2)
res = {}
for est in ('lr', 'ibp'):
    for M in (20000, 80000):
        t0 = time.time(); d = dmft.solve(T, q, M=M, seed=1, estimator=est); res[(est, M)] = d
        G1 = np.array([d['G'][t + 1, t] for t in range(S)]); cex = d['c'] / np.sqrt(2000)
        print(f"{est:3s} M={M:6d} ({time.time()-t0:.1f}s): G(t+1,t) rel err max {np.max(np.abs(G1[:180]-cex[:180])/cex[:180]):.3f}; "
              + " ".join(f"su[{t}]={d['su'][t]:.3f}" for t in (25, 50, 100, 150, 200)))
print("sim  B=512:" + " ".join(f" su[{t}]={r['su'][t]:.3f}" for t in (25, 50, 100, 150, 200)))
