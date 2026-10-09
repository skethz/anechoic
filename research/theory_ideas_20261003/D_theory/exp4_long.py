import os
for k in ("OMP_NUM_THREADS","OPENBLAS_NUM_THREADS","VECLIB_MAXIMUM_THREADS"): os.environ[k]="1"
import numpy as np, time
from sca_onsager import make_J, run
N=400; R=64
inst=[make_J(N, seed=s) for s in (1,2,3)]
best=np.array([-0.74875,-0.74125,-0.746])
t0=time.time()
for sched,(T0,S) in {"short":(30,560),"long":(40,1560)}.items():
    for lam in [0.0,0.5,0.7,0.8,0.9,1.0]:
        efs=[];trs=[]
        for i,J in enumerate(inst):
            o=run(J, lam=lam, R=R, seed=200+i, T0=T0, S=S, record=True); efs.append(o["E_final"]/N**1.5); trs.append(o["trace"])
        ef=np.stack(efs); best=np.minimum(best, ef.min(1))
        b=np.mean([tr["lam_eff"]*tr["x"] for tr in trs],0); x=np.mean([tr["x"] for tr in trs],0)
        fb=np.mean([tr["fb"] for tr in trs],0); fl=np.mean([tr["flips"].sum() for tr in trs])
        Ts=o["Ts"]/np.sqrt(N/2000)
        cross=np.argmax(b>1.0) if (b>1).any() else -1
        mid=slice(S//10, 9*S//10)
        print(f"{sched:5s} lam={lam:.1f} mean e={ef.mean():.4f} p(0.2%)={np.mean(ef<=(best*(1-0.002))[:,None]):.3f} "
              f"flips={fl:.0f} <x>mid={x[mid].mean():.2f} max b={b.max():.2f} b>1 at step {cross} (T2000={Ts[cross] if cross>=0 else float('nan'):.1f}) "
              f"flipback(mid)={fb[mid].mean():.3f}  [{time.time()-t0:.0f}s]", flush=True)
