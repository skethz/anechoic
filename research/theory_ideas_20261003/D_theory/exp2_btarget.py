# Test: is the failure controlled by b = lam * n_lin/(4T^2)?  Compare fixed-lam vs fixed-b controller.
import os
for k in ("OMP_NUM_THREADS","OPENBLAS_NUM_THREADS","VECLIB_MAXIMUM_THREADS"): os.environ[k]="1"
import numpy as np, time
from sca_onsager import make_J, run
N=400; R=96
inst=[make_J(N, seed=s) for s in (1,2,3)]
t0=time.time()
rows=[]
def summarize(tag, outs):
    ef=np.concatenate([o["E_final"]/N**1.5 for o in outs])
    bmax=np.max([np.max(o["trace"]["lam_eff"]*o["trace"]["x"]) for o in outs])
    fl=np.mean([o["trace"]["flips"].sum() for o in outs])
    print(f"{tag:22s} e_final mean={ef.mean():.4f}  median={np.median(ef):.4f}  max(b)={bmax:.2f}  flips={fl:.0f}  [{time.time()-t0:.0f}s]", flush=True)
    return ef
res={}
for lam in [0.0,0.5,0.7,0.8,0.85,0.9,1.0]:
    outs=[run(J, lam=lam, R=R, seed=100+i, record=True) for i,J in enumerate(inst)]
    res[("lam",lam)]=summarize(f"fixed lam={lam}", outs)
for b in [0.2,0.35,0.5,0.65,0.8,0.9,1.0,1.1,1.25]:
    outs=[run(J, R=R, seed=100+i, record=True, b_target=b) for i,J in enumerate(inst)]
    res[("b",b)]=summarize(f"fixed b={b}", outs)
np.save("exp2_res.npy", res, allow_pickle=True)
