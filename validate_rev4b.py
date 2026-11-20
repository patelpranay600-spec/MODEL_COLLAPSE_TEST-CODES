import numpy as np
from scipy.stats import norm
from scipy.ndimage import gaussian_filter1d
mus=[-3,2,8]; sds=[1,0.8,0.6]; pis=[0.70,0.27,0.03]; b12,b23=-0.10488,5.58249
def leak_formula(lam,h,sig,T):
    a=1-lam; out=np.zeros((3,3))
    edges=[-np.inf,b12,b23,np.inf]
    for k in range(3):
        for j in range(T+1):
            wj = lam*a**j if j<T else a**T
            sj2=h*h+j*(h*h+sig*sig)
            s=np.sqrt(sds[k]**2+sj2)
            for r in range(3):
                out[k,r]+=pis[k]*wj*(norm.cdf((edges[r+1]-mus[k])/s)-norm.cdf((edges[r]-mus[k])/s))
    return out
dx=0.01; x=np.arange(-40,55,dx)
comp=[pis[k]*norm.pdf(x,mus[k],sds[k]) for k in range(3)]
reg=[x<b12,(x>=b12)&(x<b23),x>=b23]
def grid(lam,h,T):
    a=1-lam; q=[gaussian_filter1d(c,h/dx,mode='constant') for c in comp]; R=list(q)
    for t in range(T): q=[gaussian_filter1d(a*q[k],h/dx,mode='constant')+lam*R[k] for k in range(3)]
    return np.array([[q[k][reg[r]].sum()*dx for r in range(3)] for k in range(3)])
for lam,h,T in [(0.05,0.6,59),(0.02,0.3,149),(0.05,0.6,800)]:
    F=leak_formula(lam,h,0.0,T); G=grid(lam,h,T)
    print(f"lam={lam} h={h} T={T}: max|formula-grid| = {np.abs(F-G).max():.2e}")
    print(" leakage into R3 by lineage [k=1,2,3]:",np.round(F[:,2],5)," share of own lineage:",np.round(F[:,2]/np.array(pis),4))
    print(" own-region retention A[3,3]/omega_3 =",round(F[2,2]/pis[2],4), " rare-region mass",round(F[:,2].sum(),5),"inflow",round(F[:2,2].sum(),5))
# one-step probability from component centre, for the record
for h,sig in [(0.6,0.0),(0.6,0.5)]:
    print(f"one-step P(enter R3) from x=2: h={h},sigma={sig}:", norm.sf((b23-2)/np.sqrt(h*h+sig*sig)))
