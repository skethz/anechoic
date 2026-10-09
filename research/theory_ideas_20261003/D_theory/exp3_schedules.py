# Is the lam>0 gain just an effective-temperature (schedule) shift?  Compare lam=0 with
# re-tuned schedules against lam=0.5 / 0.7 on the standard schedule; also memory kernel.
import os
for k in ("OMP_NUM_THREADS","OPENBLAS_NUM_THREADS","VECLIB_MAXIMUM_THREADS"): os.environ[k]="1"
import numpy as np, time
from sca_onsager import make_J, run
N=400; R=96
inst=[make_J(N, seed=s) for s in (1,2,3)]
best=np.array([-0.74875,-0.74125,-0.746])   # from exp2 (min over all runs)
t0=time.time()
def evaluate(tag, **kw):
    efs=[]
    for i,J in enumerate(inst):
        o=run(J, R=R, seed=100+i, **kw); efs.append(o["E_final"]/N**1.5)
    ef=np.stack(efs)
    p2=np.mean(ef <= (best*(1-0.002))[:,None]); p4=np.mean(ef <= (best*(1-0.004))[:,None])
    print(f"{tag:34s} mean e={ef.mean():.4f}  p(0.2%)={p2:.3f}  p(0.4%)={p4:.3f}  [{time.time()-t0:.0f}s]", flush=True)
    return ef
out={}
out["lam0 30-5"]=evaluate("lam=0   T 30->5", lam=0.0)
out["lam0.5 30-5"]=evaluate("lam=0.5 T 30->5", lam=0.5)
out["lam0.7 30-5"]=evaluate("lam=0.7 T 30->5", lam=0.7)
for T0,T1 in [(25,5),(35,5),(40,5),(30,4),(30,6),(36,6)]:
    out[f"lam0 {T0}-{T1}"]=evaluate(f"lam=0   T {T0}->{T1}", lam=0.0, T0=T0, T1=T1)
for q in [2.0, 3.0, 6.0]:
    out[f"lam0 q{q}"]=evaluate(f"lam=0   q={q}", lam=0.0, q=q)
for lam in [0.3,0.5,0.7]:
    out[f"mem lam{lam}"]=evaluate(f"memory-kernel lam={lam}", lam=lam, gamma_mode="qchi")
np.save("exp3_res.npy", out, allow_pickle=True)
