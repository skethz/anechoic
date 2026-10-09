# Local-field ("stay margin") distribution near the band during the anneal: does count/(4T^2) collapse
# when plotted vs y = (s h)/(2T)?  (linear pseudogap smeared on scale T  =>  n_lin ~ x0 * 4T^2)
import os
for k in ("OMP_NUM_THREADS","OPENBLAS_NUM_THREADS","VECLIB_MAXIMUM_THREADS"): os.environ[k]="1"
import numpy as np
from sca_onsager import make_J
N=400; R=64; J=make_J(N,1); sc=np.sqrt(N/2000); q=4*sc
S=560; Ts=30*sc*(5/30)**(np.arange(S)/(S-1))
bins=np.arange(-3,4.01,0.5)
for lam in (0.0,0.5):
    rng=np.random.default_rng(7); s=rng.choice([-1.,1.],size=(N,R)); sp=s.copy(); c=np.zeros(R)
    snaps={}
    for t in range(S):
        T=Ts[t]; h=J@s
        z=s*(h-lam*c*sp)+q
        if t in (100,200,300,400):
            y=(s*h)/(2*T)   # uncorrected margin in band units
            H,_=np.histogram(y,bins=bins); snaps[t]=(T/sc,H/R/(4*T*T))
        P=np.clip(z/(4*T)+0.5,0,1); fl=P<rng.random((N,R))
        c=(np.abs(z)<2*T).sum(0)/(2*T); sp=s; s=np.where(fl,-s,s)
    print(f"lam={lam}: counts per replica / (4T^2) in bins of y=s*h/(2T) (bin width 0.5)")
    print("   y-bin:   "+" ".join(f"{b:5.1f}" for b in bins[:-1]))
    for t,(T2000,H) in snaps.items():
        print(f"   T2000={T2000:5.1f} "+" ".join(f"{v:5.3f}" for v in H))
