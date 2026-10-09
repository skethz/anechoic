import os
for k in ("OMP_NUM_THREADS","OPENBLAS_NUM_THREADS","VECLIB_MAXIMUM_THREADS"): os.environ[k]="1"
import numpy as np, time, sys
from sca_onsager import make_J, run
N=400; J=make_J(N, seed=1)
lams=[-0.5,0.0,0.25,0.5,0.6,0.7,0.8,0.9,1.0]
res={}
t0=time.time()
for lam in lams:
    r=run(J, lam=lam, R=64, seed=11, record=True)
    res[lam]=r
    tr=r["trace"]
    e=r["E_final"]/N**1.5; em=r["E_min"]/N**1.5
    print(f"lam={lam:5.2f}  e_final mean={e.mean():.4f} best={e.min():.4f}  e_min mean={em.mean():.4f}  flips_tot={tr['flips'].sum():.0f}  max(lam*x)={(lam*tr['x']).max():.2f}  t={time.time()-t0:.1f}s", flush=True)
np.save("exp1_res.npy", {k:{"Ef":v["E_final"],"Em":v["E_min"],**{kk:vv for kk,vv in v["trace"].items()}} for k,v in res.items()}, allow_pickle=True)
