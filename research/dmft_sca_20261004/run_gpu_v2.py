import json, sys, time
import numpy as np
import dmft_torch as D
def sched(t0, S): return t0 * (5.0 / t0) ** (np.arange(S) / max(1, S - 1))
def log(msg): print(f"{time.strftime('%H:%M:%S')} {msg}", flush=True)
M = int(sys.argv[1]); itw = int(sys.argv[2]); cases = sys.argv[3].split(','); seed = int(sys.argv[4]) if len(sys.argv) > 4 else 21
spec = {'S200_plain': (200, 30.0, 0.0), 'S560_lam0': (560, 30.0, 0.0), 'S560_lam0.7': (560, 30.0, 0.7)}
out = {}
for name in cases:
    S, t0, lam = spec[name]; t1 = time.time()
    d = D.solve_causal2(sched(t0, S), np.full(S, 4.0), lam=(np.full(S, lam) if lam else None), M=M, it_window=itw, seed=seed)
    pts = (25, 50, 100, 150, 200) if S == 200 else (140, 280, 420, 560)
    out[name] = dict(su=d['su'].tolist(), flips=d['flips'].tolist(), plin=d['plin'].tolist(),
                     C_lag1=[float(d['C'][t, t - 1]) for t in range(1, S + 1)],
                     G_lags=[[float(d['G'][t, t - k]) if t - k >= 0 else 0.0 for k in range(1, 9)] for t in range(S + 1)])
    log(f"v2 {name} M{M} it{itw} s{seed}: su " + " ".join(f"{d['su'][t]:.4f}" for t in pts) + f" | flips/run {sum(d['flips'])*2000:.0f} | {time.time()-t1:.0f}s")
    json.dump(out, open(f'results/gpu_v2_M{M}_it{itw}_s{seed}_{"_".join(cases)}.json', 'w'))
