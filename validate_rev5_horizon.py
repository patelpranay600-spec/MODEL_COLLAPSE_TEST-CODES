# R5-1b: sustained 5%-of-gap stabilisation time vs T5 across (lambda,h); deterministic uncorrected target-averaged operator
import numpy as np
from scipy.stats import norm
from scipy.ndimage import gaussian_filter1d
mus=[-3,2,8]; sds=[1,0.8,0.6]; pis=[0.70,0.27,0.03]; b12,b23=-0.10488,5.58249
dx=0.02; x=np.arange(-40,60,dx)
comp=[pis[k]*norm.pdf(x,mus[k],sds[k]) for k in range(3)]; P0=sum(comp)
reg2=x>=b23; tail=(x<=-5.5)|(x>=9.5); F0=np.cumsum(P0)*dx
tau0=P0[tail].sum()*dx; mR0=0.03
def ep(q):
    Q=sum(q); mR=Q[reg2].sum()*dx; tau=Q[tail].sum()*dx
    return [abs(mR-mR0)/mR0, abs(np.log(tau/tau0)), np.abs(np.cumsum(Q)*dx-F0).sum()*dx, q[2][reg2].sum()*dx, (q[0][reg2].sum()+q[1][reg2].sum())*dx]
names=['e_R','e_tau','W1','A33','inflow']
def stab(lam,h):
    a=1-lam; T5=int(np.ceil(np.log(0.05)/np.log(a))); Tmax=int(np.ceil(np.log(1e-9)/np.log(a)))
    q=[gaussian_filter1d(c,h/dx,mode='constant') for c in comp]; R=list(q)
    tr=[ep(q)]
    for t in range(Tmax):
        q=[gaussian_filter1d(a*q[k],h/dx,mode='constant')+lam*R[k] for k in range(3)]; tr.append(ep(q))
    tr=np.array(tr); out=[]
    for j in range(5):
        v=tr[:,j]; E0,Es=v[0],v[-1]; ok=np.abs(v-Es)<=0.05*abs(Es-E0)
        t=next((t for t in range(len(ok)) if ok[t:].all()),None); out.append(t)
    return T5,out
print("t_stab = first t with |E_t' - E*| <= 5% of |E* - E_0| for all t'>=t   (grid dx=0.02; E* at residual 1e-9)")
print("lam    h    T5 | "+" ".join(f"{n:>7s}" for n in names)+" | max over (e_R,e_tau,W1) / T5")
for lam in [0.01,0.02,0.05,0.10,0.20]:
    for h in [0.3,0.6]:
        T5,o=stab(lam,h)
        print(f"{lam:4.2f} {h:4.1f} {T5:4d} | "+" ".join(f"{t:7d}" for t in o)+f" | {max(o[:3])/T5:.2f}",flush=True)
