# Test the band-shift self-consistency prediction  x = x0 (1 + lam x)^2  ->  lam_c = 1/(4 x0)
# for several inertia values q and for Gaussian SK couplings.
import os
for k in ("OMP_NUM_THREADS","OPENBLAS_NUM_THREADS","VECLIB_MAXIMUM_THREADS"): os.environ[k]="1"
import numpy as np, time
from sca_onsager import make_J, run
N=400; R=48; t0=time.time()
def xmid(tr,S): return tr["x"][S//10:8*S//10].mean()
def xpred(lam,x0):
    a=lam*x0
    if 4*a>1: return np.nan
    return (1-2*a-np.sqrt(1-4*a))/(2*x0*lam*lam) if lam>0 else x0
for kind in ("pm1","gauss"):
    inst=[make_J(N, seed=s, kind=kind) for s in (11,12)]
    for q in ([4.0,8.0,12.0] if kind=="pm1" else [4.0]):
        rows=[]
        for lam in [0.0,0.3,0.5,0.6,0.7,0.8,0.9,1.0,1.2,1.4]:
            outs=[run(J, lam=lam, R=R, seed=300+i, q=q, record=True) for i,J in enumerate(inst)]
            e=np.mean([o["E_final"].mean() for o in outs])/N**1.5
            xm=np.mean([xmid(o["trace"],560) for o in outs])
            bmax=np.max([np.max(o["trace"]["lam_eff"]*o["trace"]["x"]) for o in outs])
            fl=np.mean([o["trace"]["flips"].sum() for o in outs])
            rows.append((lam,e,xm,bmax,fl))
        x0=rows[0][2]; fl0=rows[0][4]
        print(f"== {kind} q={q}: x0={x0:.3f} -> predicted lam_c=1/(4x0)={1/(4*x0):.2f}   [{time.time()-t0:.0f}s]")
        for lam,e,xm,bmax,fl in rows:
            print(f"   lam={lam:.1f} e={e:.4f} x_mid={xm:.3f} (pred {xpred(lam,x0):.3f}) max b={bmax:.2f} flips/flips0={fl/fl0:.2f} (pred (1+b)^2={(1+lam*xpred(lam,x0))**2 if lam*x0<=0.25 else float('nan'):.2f})", flush=True)
