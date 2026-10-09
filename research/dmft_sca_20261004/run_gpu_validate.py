"""GPU DMFT validation: same three cases as the CPU solver (S=200 plain; S=560 lambda 0 and 0.7) at large M."""
import json, sys, time
import numpy as np
import dmft_torch as D

T_FIN = 5.0


def sched(t0, S):
    return t0 * (T_FIN / t0) ** (np.arange(S) / max(1, S - 1))


def log(msg):
    print(f"{time.strftime('%H:%M:%S')} {msg}", flush=True)


M = int(sys.argv[1]) if len(sys.argv) > 1 else 1_000_000
itw = int(sys.argv[2]) if len(sys.argv) > 2 else 32
out = {}
cases = [('S200_plain', 200, 30.0, 0.0), ('S560_lam0', 560, 30.0, 0.0), ('S560_lam0.7', 560, 30.0, 0.7)]
for name, S, t0, lam in cases:
    t1 = time.time()
    d = D.solve_causal(sched(t0, S), np.full(S, 4.0), lam=(np.full(S, lam) if lam else None), M=M, it_window=itw, seed=11, log=log)
    out[name] = dict(su=d['su'].tolist(), flips_per_run=float(d['flips'].sum() * 2000), plin=d['plin'].tolist(), c=d['c'].tolist(),
                     G_lags=[[float(d['G'][t, t - k]) if t - k >= 0 else 0.0 for k in range(1, 9)] for t in range(S + 1)],
                     sec=time.time() - t1, M=M, it_window=itw)
    pts = (25, 50, 100, 150, 200) if S == 200 else (140, 280, 420, 560)
    log(f"{name}: su " + " ".join(f"{d['su'][t]:.4f}" for t in pts) + f" | flips/run {out[name]['flips_per_run']:.0f} | {time.time()-t1:.0f}s")
    json.dump(out, open(f'results/gpu_validate_M{M}_it{itw}.json', 'w'))
