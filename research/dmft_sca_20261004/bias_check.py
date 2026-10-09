import sys, time, numpy as np
sys.path.insert(0, str(__import__('pathlib').Path(__file__).resolve().parent / '../theory_ideas_20261003'))
import sweeps_abc as m
import dmft
S = 200; T = m.sched(30.0, S); q = np.full(S, 4.0)
print("ref (N=4000 sims)            0.7930 1.1680 1.4236 1.4698 1.4778")
for M, itw in ((30000, 16), (120000, 16), (30000, 64)):
    t0 = time.time(); d = dmft.solve_causal(T, q, M=M, it_window=itw, seed=5)
    print(f"causal M={M:6d} it_window={itw:2d} ({time.time()-t0:4.0f}s): " + " ".join(f"{d['su'][t]:.4f}" for t in (25, 50, 100, 150, 200)), flush=True)
