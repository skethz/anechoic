import sys, time, json, numpy as np
sys.path.insert(0, str(__import__('pathlib').Path(__file__).resolve().parent / '../theory_ideas_20261003'))
import sweeps_abc as m
import dmft
lam_sim = json.load(open('lambda_scan.json'))
# S=200 plain reference (random-instance mean, N=4000): su[25,50,100,150,200] = 0.793 1.168 1.4236 1.4698 1.4778
S = 200; T = m.sched(30.0, S); q = np.full(S, 4.0)
for seed in (1, 2):
    t0 = time.time(); d = dmft.solve_causal(T, q, M=30000, seed=seed)
    print(f"S200 plain s{seed} ({time.time()-t0:.0f}s): " + " ".join(f"{d['su'][t]:.4f}" for t in (25, 50, 100, 150, 200))
          + "  [ref 0.793 1.168 1.4236 1.4698 1.4778]", flush=True)
S = 560; T = m.sched(30.0, S); q = np.full(S, 4.0)
for lam in (0.0, 0.7):
    ref = lam_sim[str(lam)]['sim_su']
    for seed in (1, 2):
        t0 = time.time(); d = dmft.solve_causal(T, q, lam=np.full(S, lam), M=30000, seed=seed)
        print(f"S560 lambda {lam} s{seed} ({time.time()-t0:.0f}s): su[140,280,420,560] " + " ".join(f"{d['su'][t]:.4f}" for t in (140, 280, 420, 560))
              + "  [sim " + " ".join(f"{ref[t]:.4f}" for t in (140, 280, 420, 560)) + "]", flush=True)
