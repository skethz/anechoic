"""Does DMFT reproduce the Onsager-SCA threshold? Short A2 schedule: S=560, T 30->5, q=4, constant lambda (A2 P3: p(0.8)=0.612, p(0.9)=0.235, p(1.0)=0)."""
import sys, time, json, numpy as np
sys.path.insert(0, str(__import__('pathlib').Path(__file__).resolve().parent / '../theory_ideas_20261003'))
import sweeps_abc as m
import dmft
J, sumw = m.load()
S = 560; T = m.sched(30.0, S); q = np.full(S, 4.0)
res = {}
for lam in (0.0, 0.5, 0.7, 0.8, 0.9, 1.0):
    L = np.full(S, lam)
    t0 = time.time()
    d = dmft.solve_iter(T, q, lam=L, M=30000, iters=60, alpha=0.25, avg_last=20, seed=7)
    r = dmft.simulate(J, T, q, lam=L, B=256, seed=8)
    res[lam] = dict(dmft_su=d['su'].tolist(), sim_su=r['su'].tolist(), dmft_flips=float(d['flips'].sum() * 2000),
                    sim_flips=float(r['flips'].sum() * 2000), sim_cut=float(r['cuts'].mean()), sim_p=float((r['cuts'] >= 33000).mean()))
    print(f"lambda {lam:.1f} ({time.time()-t0:.0f}s): final su dmft {d['su'][-1]:.4f} sim {r['su'][-1]:.4f} | "
          f"cut dmft {dmft.cut_from_su(d['su'][-1]):.0f} sim {r['cuts'].mean():.0f} (p_sim {res[lam]['sim_p']:.3f}) | "
          f"flips/run dmft {res[lam]['dmft_flips']:.0f} sim {res[lam]['sim_flips']:.0f}", flush=True)
json.dump(res, open('lambda_scan.json', 'w'))
