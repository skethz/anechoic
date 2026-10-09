import sys, time, json, numpy as np
sys.path.insert(0, str(__import__('pathlib').Path(__file__).resolve().parent / '../theory_ideas_20261003'))
import sweeps_abc as m
import dmft
S = 200; T = m.sched(30.0, S); q = np.full(S, 4.0)
pts = (25, 50, 100, 150, 200); res = {}
for seed in (1, 2):
    t0 = time.time(); d = dmft.solve_iter(T, q, M=50000, iters=120, alpha=0.2, avg_last=40, seed=seed)
    tr = [x['su_final'] for x in d['hist']]
    print(f"s{seed} ({time.time()-t0:.0f}s): " + " ".join(f"{d['su'][t]:.4f}" for t in pts)
          + " | su_final iters 40..120 step 10: " + " ".join(f"{np.mean(tr[i:i+10]):.3f}" for i in range(40, 120, 10)), flush=True)
    res[seed] = {k: (v.tolist() if hasattr(v, 'tolist') else v) for k, v in d.items() if k in ('su', 'flips', 'plin', 'c')}
json.dump(res, open('conv_long.json', 'w'))
