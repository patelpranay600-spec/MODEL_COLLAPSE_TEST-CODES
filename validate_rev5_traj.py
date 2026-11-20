# R5-1: endpoint trajectories on the deterministic grid (target-averaged operator, sigma=0)
import numpy as np
from scipy.stats import norm
from scipy.ndimage import gaussian_filter1d
mus=[-3,2,8]; sds=[1,0.8,0.6]; pis=[0.70,0.27,0.03]; b12,b23=-0.10488,5.58249
dx=0.01; x=np.arange(-40,60,dx)
comp=[pis[k]*norm.pdf(x,mus[k],sds[k]) for k in range(3)]; P0=sum(comp)
reg=[x<b12,(x>=b12)&(x<b23),x>=b23]; tail=(x<=-5.5)|(x>=9.5)
F0=np.cumsum(P0)*dx
tau0=P0[tail].sum()*dx; mR0=0.03
def endpoints(q):
    Q=sum(q); mR=Q[reg[2]].sum()*dx; tau=Q[tail].sum()*dx
    W1=np.abs(np.cumsum(Q)*dx-F0).sum()*dx
    A33=q[2][reg[2]].sum()*dx; inflow=q[0][reg[2]].sum()*dx+q[1][reg[2]].sum()*dx
    return dict(eR=abs(mR-mR0)/mR0, epsR=(mR-mR0)/mR0, etau=abs(np.log(tau/tau0)), W1=W1, A33=A33, inflow=inflow)
def traj(lam,h,Tmax):
    a=1-lam; q=[gaussian_filter1d(c,h/dx,mode='constant') for c in comp]; R=list(q)
    out=[endpoints(q)]
    for t in range(1,Tmax+1):
        q=[gaussian_filter1d(a*q[k],h/dx,mode='constant')+lam*R[k] for k in range(3)]
        out.append(endpoints(q))
    return out
def sustained(vals,Estar,E0,rel=0.05):
    thr=rel*abs(Estar-E0)
    ok=np.abs(np.array(vals)-Estar)<=thr
    # first t such that ok for all t'>=t
    for t in range(len(ok)):
        if ok[t:].all(): return t
    return None
for lam,h,Tmax in [(0.05,0.6,800),(0.02,0.3,1500)]:
    T5=int(np.ceil(np.log(0.05)/np.log(1-lam)))
    tr=traj(lam,h,Tmax)
    print(f"\n=== lam={lam} h={h}  T5={T5}  (E* taken at t={Tmax}) ===")
    keys=['eR','epsR','etau','W1','A33','inflow']
    print("t".rjust(5)," ".join(k.rjust(9) for k in keys))
    for t in sorted(set([0,5,10,20,T5,2*T5,3*T5,Tmax])):
        if t<=Tmax: print(str(t).rjust(5)," ".join(f"{tr[t][k]:9.5f}" for k in keys))
    print("fraction of gap (E0->E*) closed at t=T5 / 2*T5, and sustained 5%-crossing time:")
    for k in keys:
        v=[d[k] for d in tr]; E0,Es=v[0],v[-1]
        gap=Es-E0
        f=lambda t: (v[t]-E0)/gap if abs(gap)>1e-12 else float('nan')
        print(f"  {k:7s} E0={E0:9.5f} E*={Es:9.5f}  closed@T5={f(T5):6.3f} @2T5={f(2*T5):6.3f}  t_stab={sustained(v,Es,E0)}  value@T5-4..T5 mean={np.mean(v[T5-4:T5+1]):9.5f}")
