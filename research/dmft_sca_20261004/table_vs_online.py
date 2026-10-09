"""Theory prediction 5: c(t) = n_lin/(2T) is self-averaging, so a precomputed table should equal the online popcount.
O1 configuration; table = mean online c(t) from a 256-run pilot (different seeds); 2048 paired runs per variant."""
import json, sys, time
import numpy as np
sys.path.insert(0, str(__import__('pathlib').Path(__file__).resolve().parent / '../theory_ideas_20261003'))
import sweeps_abc as m
import multilag as ML
J, sumw = m.load(); N = len(J)
S, t0, q, lam = 960, 12.0, 8.0, 1.05
T = m.sched(t0, S); L1 = ML.lam_sched(lam, True, S)


def run(B, rng, table=None, record=False):
    s = rng.choice(np.array([-1.0, 1.0], np.float32), size=(B, N)); h = s @ J; sp = None; c_prev = None
    flips = np.zeros(B); crec = np.zeros(S)
    for t in range(S):
        field = h.copy()
        if sp is not None:
            c = table[t] if table is not None else c_prev
            field -= L1[t] * c * sp
        z = s * field + q
        flip = np.clip(z / (4 * T[t]) + 0.5, 0, 1) < rng.random((B, N), dtype=np.float32)
        c_prev = ((z > -2 * T[t]) & (z < 2 * T[t])).sum(1, keepdims=True) / (2 * T[t])
        if record and t + 1 < S: crec[t + 1] = float(c_prev.mean())
        d = np.where(flip, -2 * s, 0).astype(np.float32)
        sp = s; s = s + d; h = h + d @ J; flips += flip.sum(1)
    return m.cut(J, sumw, s), flips, crec


_, _, table = run(256, np.random.default_rng(5), record=True)
cv = np.array([run(1, np.random.default_rng(100 + i), record=True)[2] for i in range(4)])  # per-run spread of c (B=1)
print("c(t) table at t=100,300,600,900:", " ".join(f"{table[t]:.2f}" for t in (100, 300, 600, 900)),
      "| single-run c spread (sd/mean) at same t:", " ".join(f"{cv[:, t].std()/max(cv[:, t].mean(),1e-9):.3f}" for t in (100, 300, 600, 900)))
res = {}
for name, tab in (('online', None), ('table', table)):
    t1 = time.time(); cut, fl, _ = run(2048, np.random.default_rng(4242), table=tab)
    p = float((cut >= 33000).mean()); res[name] = dict(p=p, mean_cut=float(cut.mean()), flips=float(fl.mean()))
    print(f"{name:6s}: p={p:.3f} mean_cut={cut.mean():.1f} flips={fl.mean():.0f} ({time.time()-t1:.0f}s)", flush=True)
json.dump(dict(table=table.tolist(), res=res), open('table_vs_online.json', 'w'))
