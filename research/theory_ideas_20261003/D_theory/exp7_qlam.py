# Joint (q, lam) test + schedule rules, with the same instances/thresholds as exp5.
import os
for k in ("OMP_NUM_THREADS","OPENBLAS_NUM_THREADS","VECLIB_MAXIMUM_THREADS"): os.environ[k]="1"
import numpy as np, time
from sca_onsager import make_J, run
N=400; R=192; S=560
inst=[make_J(N, seed=s) for s in (1,2,3,4)]
prev=np.load("exp5_res.npy",allow_pickle=True).item(); best=prev["best"].copy()
t0=time.time(); E={}
cfgs={
 "q4 lam0":dict(q=4.0,lam=0.0),
 "q4 lam0.6":dict(q=4.0,lam=0.6),
 "q8 lam0":dict(q=8.0,lam=0.0),
 "q8 lam0.7":dict(q=8.0,lam=0.7),
 "q8 lam0.9":dict(q=8.0,lam=0.9),
 "q8 lam1.0":dict(q=8.0,lam=1.0),
 "q12 lam0.9":dict(q=12.0,lam=0.9),
 "q4 lam0.7->0 (last 30%)":dict(q=4.0,lam_sched=lambda t,T: 0.7 if t<0.7*S else 0.7*(S-t)/(0.3*S)),
 "q4 capped b=0.3, lam<=0.7":"bcap",
}
for name,kw in cfgs.items():
    if kw=="bcap":
        # lam(t)=min(0.7, 0.3/x0_est) using x from the lam=0 pilot is equivalent to a constant; emulate online:
        kw=dict(q=4.0,lam_sched=None,b_target=None)
        efs=[]
        for i,J in enumerate(inst):
            # pilot to get x0 (median over mid-anneal), then fixed lam = b*/(x0(1+b*)^2)
            pil=run(J, lam=0.0, R=16, seed=999+i, record=True)
            x0=np.median(pil["trace"]["x"][S//10:8*S//10]); bstar=0.4
            lam=bstar/(x0*(1+bstar)**2)
            efs.append(run(J, lam=lam, R=R, seed=2000+7*i, q=4.0)["E_final"]/N**1.5)
            print(f"   inst {i}: x0={x0:.3f} -> lam={lam:.3f} (lam_c=1/(4x0)={1/(4*x0):.2f})")
        E[name]=np.stack(efs)
    else:
        E[name]=np.stack([run(J, R=R, seed=2000+7*i, **kw)["E_final"]/N**1.5 for i,J in enumerate(inst)])
    print(f"{name:28s} done [{time.time()-t0:.0f}s]", flush=True)
best=np.minimum(best, np.min(np.stack([v.min(1) for v in E.values()]),0))
print("best",best)
for name,ef in E.items():
    ps=[np.mean(ef<=(best*(1-d))[:,None]) for d in (0.002,0.004,0.008)]
    print(f"{name:28s} mean e={ef.mean():.4f}±{ef.std()/np.sqrt(ef.size):.4f}  p0.2%={ps[0]:.3f}  p0.4%={ps[1]:.3f}  p0.8%={ps[2]:.3f}")
np.save("exp7_res.npy",{"E":E,"best":best},allow_pickle=True)
