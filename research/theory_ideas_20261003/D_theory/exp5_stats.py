# Careful statistics: success prob vs lam at fixed thresholds (relative to best-found per instance).
import os
for k in ("OMP_NUM_THREADS","OPENBLAS_NUM_THREADS","VECLIB_MAXIMUM_THREADS"): os.environ[k]="1"
import numpy as np, time
from sca_onsager import make_J, run
N=400; R=192
inst=[make_J(N, seed=s) for s in (1,2,3,4)]
lams=[0.0,0.25,0.5,0.6,0.7,0.8]
t0=time.time(); E={}
for S in (250,560):
    for lam in lams:
        E[(S,lam)]=np.stack([run(J, lam=lam, R=R, seed=1000+7*i+S, S=S)["E_final"]/N**1.5 for i,J in enumerate(inst)])
    print("S",S,"done",time.time()-t0,flush=True)
best=np.min(np.stack([v.min(1) for v in E.values()]),0)
print("best per instance",best)
np.save("exp5_res.npy",{"E":E,"best":best},allow_pickle=True)
for S in (250,560):
    print(f"--- S={S}")
    for lam in lams:
        ef=E[(S,lam)]
        ps=[np.mean(ef<=(best*(1-d))[:,None]) for d in (0.002,0.004,0.008)]
        se=[np.sqrt(p*(1-p)/ef.size) for p in ps]
        print(f"lam={lam:.2f} mean e={ef.mean():.4f}±{ef.std()/np.sqrt(ef.size):.4f}  p0.2%={ps[0]:.3f}±{se[0]:.3f}  p0.4%={ps[1]:.3f}±{se[1]:.3f}  p0.8%={ps[2]:.3f}±{se[2]:.3f}")
