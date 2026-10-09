# Ramp-down of lam at the end of the anneal, combined with larger q.
import os
for k in ("OMP_NUM_THREADS","OPENBLAS_NUM_THREADS","VECLIB_MAXIMUM_THREADS"): os.environ[k]="1"
import numpy as np, time
from sca_onsager import make_J, run
N=400; R=192; S=560
inst=[make_J(N, seed=s) for s in (1,2,3,4)]
best=np.load("exp7_res.npy",allow_pickle=True).item()["best"]
def ramp(l0,frac=0.3):
    return lambda t,T: l0 if t<(1-frac)*S else l0*(S-t)/(frac*S)
cfgs={
 "q4 lam0.5 ramp30%":dict(q=4.0,lam_sched=ramp(0.5)),
 "q4 lam0.7 ramp50%":dict(q=4.0,lam_sched=ramp(0.7,0.5)),
 "q4 lam0.8 ramp30%":dict(q=4.0,lam_sched=ramp(0.8)),
 "q8 lam0.9 ramp30%":dict(q=8.0,lam_sched=ramp(0.9)),
 "q8 lam1.0 ramp30%":dict(q=8.0,lam_sched=ramp(1.0)),
 "q4 lam0.7 off@70% (step)":dict(q=4.0,lam_sched=lambda t,T: 0.7 if t<0.7*S else 0.0),
}
t0=time.time(); E={}
for name,kw in cfgs.items():
    E[name]=np.stack([run(J, R=R, seed=2000+7*i, **kw)["E_final"]/N**1.5 for i,J in enumerate(inst)])
    ef=E[name]; ps=[np.mean(ef<=(best*(1-d))[:,None]) for d in (0.002,0.004,0.008)]
    print(f"{name:26s} mean e={ef.mean():.4f}±{ef.std()/np.sqrt(ef.size):.4f}  p0.2%={ps[0]:.3f}  p0.4%={ps[1]:.3f}  p0.8%={ps[2]:.3f}  [{time.time()-t0:.0f}s]",flush=True)
print("best now", np.minimum(best,np.min(np.stack([v.min(1) for v in E.values()]),0)))
