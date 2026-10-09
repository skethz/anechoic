import sys, time, json, numpy as np
sys.path.insert(0, str(__import__('pathlib').Path(__file__).resolve().parent / '../theory_ideas_20261003'))
import sweeps_abc as m
import dmft
lam_sim = json.load(open('lambda_scan.json'))
S = 560; T = m.sched(30.0, S); q = np.full(S, 4.0); out = {}
for lam in (0.0, 0.7):
    ref = lam_sim[str(lam)]['sim_su']
    for seed in (1, 2):
        t0 = time.time(); d = dmft.solve_causal(T, q, lam=np.full(S, lam), M=30000, it_window=64, seed=seed)
        out[f"{lam}_{seed}"] = d['su'].tolist()
        print(f"S560 lambda {lam} s{seed} ({time.time()-t0:.0f}s): su[140,280,420,560] " + " ".join(f"{d['su'][t]:.4f}" for t in (140, 280, 420, 560))
              + "  [sim " + " ".join(f"{ref[t]:.4f}" for t in (140, 280, 420, 560)) + "]"
              + f"  flips/run {d['flips'].sum()*2000:.0f} [sim {lam_sim[str(lam)]['sim_flips']:.0f}]", flush=True)
json.dump(out, open('conv_causal64.json', 'w'))
