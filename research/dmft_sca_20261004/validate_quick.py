import sys, time, numpy as np
sys.path.insert(0, str(__import__('pathlib').Path(__file__).resolve().parent / '../theory_ideas_20261003'))
import sweeps_abc as m
import dmft
J, sumw = m.load()
S = 200; T = m.sched(30.0, S); q = np.full(S, 4.0)
t0 = time.time(); d = dmft.solve(T, q, M=20000, seed=1); t1 = time.time()
r = dmft.simulate(J, T, q, B=256, seed=2); t2 = time.time()
print(f"dmft {t1-t0:.1f}s  sim {t2-t1:.1f}s")
# LR-estimated G(t+1,t) vs exact P_lin/(2T) (rescaled)
G1 = np.array([d['G'][t + 1, t] for t in range(S)]); cex = d['c'] / np.sqrt(2000)
print("G(t+1,t) LR vs exact: max rel err %.3f (t<180)" % np.max(np.abs(G1[:180] - cex[:180]) / cex[:180]))
for t in (0, 10, 25, 50, 100, 150, 199, 200):
    line = f"t={t:3d} su dmft {d['su'][t]:.4f} sim {r['su'][t]:.4f}"
    if t < S:
        line += f" | flips {d['flips'][t]:.4f} {r['flips'][t]:.4f} | plin {d['plin'][t]:.4f} {r['plin'][t]:.4f}"
    print(line)
print("final cut: dmft %.0f  sim mean %.0f (sd %.0f)" % (dmft.cut_from_su(d['su'][S]), r['cuts'].mean(), r['cuts'].std()))
