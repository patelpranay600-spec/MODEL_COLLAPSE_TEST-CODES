# Experiment 1 (population vs finite-N decomposition at fixed h) and Experiment 2 (bandwidth rules), 1-D three-component GMM.
# Separate seed namespace (770xx/880xx/890xx). Replace-style recursion, N=2000, lam=0.05, with-replacement reference draw (n_ref=100), T=59.
import numpy as np, time, sys
from scipy.stats import norm
from scipy.ndimage import gaussian_filter1d
MU=np.array([-3.,2.,8.]); SD=np.array([1.,.8,.6]); PI=np.array([.70,.27,.03]); B23=-0.0+5.58249
N=2000; LAM=0.05; A_=1-LAM; NR=int(round(LAM*N)); NS=N-NR; T=59; REC=[20,55,56,57,58,59]
MS=[20,50,100,200,400]; NSEED=40; NSEL=20; ARCH=400
x0,x1,dx=-20.,30.,0.01; nG=int((x1-x0)/dx); edges=x0+dx*np.arange(nG+1); right=edges[1:]; cent=edges[:-1]+dx/2
mix_cdf=lambda t: sum(PI[k]*norm.cdf(t,MU[k],SD[k]) for k in range(3))
F0=mix_cdf(right); P0mass=np.diff(np.concatenate([[mix_cdf(edges[0])],F0]))
V0=0.70*(1+9)+0.27*(0.64+4)+0.03*(0.36+64)-(-1.32)**2
TAU0=mix_cdf(-5.5)+1-mix_cdf(9.5); MR0=1-mix_cdf(B23)
def draw_P0(rng,shape):
    k=rng.choice(3,size=shape,p=PI); return rng.normal(MU[k],SD[k]),k
def W1(cum): return np.abs(cum-F0).sum()*dx
def binpts(X): return np.bincount(np.clip(((X-x0)/dx).astype(int),0,nG-1),minlength=nG)/len(X)
def filt(mass,h): return gaussian_filter1d(mass,h/dx,mode='constant',truncate=8.0)
def pop_metrics(mass):
    cum=np.cumsum(mass); mR=1-np.interp(B23,right,cum); tau=np.interp(-5.5,right,cum)+1-np.interp(9.5,right,cum)
    mean=(mass*cent).sum(); var=(mass*cent**2).sum()-mean**2
    return dict(W1=W1(cum),epsR=(mR-MR0)/MR0,ltau=np.log(tau/TAU0) if tau>0 else np.nan,rv=var/V0,tau=tau)
def kde_metrics(X,h):   # exact region/tail/variance from points, W1 from grid
    mR=norm.sf((B23-X)/h).mean(); tau=(norm.cdf((-5.5-X)/h)+norm.sf((9.5-X)/h)).mean()
    cum=np.cumsum(filt(binpts(X),h))
    return dict(W1=W1(cum),epsR=(mR-MR0)/MR0,ltau=np.log(tau/TAU0) if tau>0 else np.nan,rv=(X.var()+h*h)/V0,tau=tau)
def raw_metrics(X):
    tau=((X<=-5.5)|(X>=9.5)).mean(); mR=(X>=B23).mean(); cum=np.cumsum(binpts(X))
    return dict(W1=W1(cum),epsR=(mR-MR0)/MR0,ltau=np.log(tau/TAU0) if tau>0 else np.nan,rv=X.var(ddof=1)/V0,tau=tau)
def pop_traj(Rmass,h):
    q=filt(A_*P0mass+LAM*Rmass,h); Rh=filt(Rmass,h); out={0:pop_metrics(q)}
    for t in range(1,T+1):
        q=filt(A_*q,h)+LAM*Rh
        if t in REC: out[t]=pop_metrics(q)
    return out
def run_chains(ref_x,ref_k,h,seed):
    C,m=ref_x.shape; rng=np.random.default_rng(seed)
    sx,sk=draw_P0(rng,(C,NS)); ri=rng.integers(0,m,(C,NR))
    X=np.concatenate([sx,np.take_along_axis(ref_x,ri,1)],1); K=np.concatenate([sk,np.take_along_axis(ref_k,ri,1)],1); rec={}
    for t in range(1,T+1):
        idx=rng.integers(0,N,(C,NS)); y=np.take_along_axis(X,idx,1)+h[:,None]*rng.standard_normal((C,NS)); yk=np.take_along_axis(K,idx,1)
        ri=rng.integers(0,m,(C,NR)); X=np.concatenate([y,np.take_along_axis(ref_x,ri,1)],1); K=np.concatenate([yk,np.take_along_axis(ref_k,ri,1)],1)
        if t in REC: rec[t]=(X.copy(),K.copy())
    return rec
def chain_metrics(rec,h):
    C=len(h); out={}
    for t,(X,K) in rec.items():
        out[t]=[dict(**kde_metrics(X[c],h[c]),omega=(K[c]==2).mean(),surv=float((K[c]==2).any()),emp0=float(((X[c]<=-5.5)|(X[c]>=9.5)).mean()==0)) for c in range(C)]
    return out
KEYS=['W1','epsR','ltau','rv']
def agg(per_t,key,spec):   # per_t: {t: list of dict}; returns array over chains
    if spec=='t20': return np.array([d[key] for d in per_t[20]])
    if spec=='t59': return np.array([d[key] for d in per_t[59]])
    return np.mean([[d[key] for d in per_t[t]] for t in [55,56,57,58,59]],0)
def aggpop(trajs,key,spec):
    if spec=='t20': return np.array([tr[20][key] for tr in trajs])
    if spec=='t59': return np.array([tr[59][key] for tr in trajs])
    return np.array([np.mean([tr[t][key] for t in [55,56,57,58,59]]) for tr in trajs])
def ms(a): a=np.asarray(a,float); a=a[~np.isnan(a)]; return a.mean(), (a.std(ddof=1)/np.sqrt(len(a)) if len(a)>1 else np.nan)
def fm(a): m,s=ms(a); return f"{m:8.4f}"
def fms(a): m,s=ms(a); return f"{m:+8.4f}±{s:.4f}"
# archives (labelled P0 draws), nested reference sets = prefixes
arch_x,arch_k=draw_P0(np.random.default_rng(77001),(NSEED,ARCH)); sel_x,sel_k=draw_P0(np.random.default_rng(77002),(NSEL,ARCH))

# ======== additions for Phase 1 ========
def run2(ref_x,ref_k,h,seed,init='P0',rec_times=(59,),track=False):
    C,m=ref_x.shape; rng=np.random.default_rng(seed)
    if init=='P0':
        sx,sk=draw_P0(rng,(C,NS)); ri=rng.integers(0,m,(C,NR))
        X=np.concatenate([sx,np.take_along_axis(ref_x,ri,1)],1); K=np.concatenate([sk,np.take_along_axis(ref_k,ri,1)],1)
    else:
        ri=rng.integers(0,m,(C,N)); X=np.take_along_axis(ref_x,ri,1); K=np.take_along_axis(ref_k,ri,1)
    om=[(K==2).mean(1)] if track else None; rec={}
    for t in range(1,T+1):
        idx=rng.integers(0,N,(C,NS)); y=np.take_along_axis(X,idx,1)+h[:,None]*rng.standard_normal((C,NS)); yk=np.take_along_axis(K,idx,1)
        ri=rng.integers(0,m,(C,NR)); X=np.concatenate([y,np.take_along_axis(ref_x,ri,1)],1); K=np.concatenate([yk,np.take_along_axis(ref_k,ri,1)],1)
        if track: om.append((K==2).mean(1))
        if t in rec_times: rec[t]=(X.copy(),K.copy())
    return rec,(np.array(om).T if track else None)
def Kcomp(t,h,init):
    if init=='ref': return [(LAM*A_**j,(j+1)*h*h) for j in range(t)]+[(A_**t,(t+1)*h*h)], None
    return [(LAM*A_**j,(j+1)*h*h) for j in range(t+1)], (A_**(t+1),(t+1)*h*h)
def Rmass_of(r): return np.bincount(np.clip(((r-x0)/dx).astype(int),0,nG-1),minlength=nG)/len(r)
def B_grid(Rmass,t,h,init):
    comps,p0=Kcomp(t,h,init); q=sum(w*filt(Rmass,np.sqrt(v)) for w,v in comps)
    return q+(p0[0]*filt(P0mass,np.sqrt(p0[1])) if p0 else 0)
def A_grid(Rmass,t,h,init):
    q=filt(Rmass,h) if init=='ref' else filt(A_*P0mass+LAM*Rmass,h); Rh=filt(Rmass,h)
    for _ in range(t): q=filt(A_*q,h)+LAM*Rh
    return q
def B_exact_cdf(r,t,h,init):
    comps,p0=Kcomp(t,h,init); F=np.zeros(nG)
    for w,v in comps: F+=w*norm.cdf((right[None,:]-r[:,None])/np.sqrt(v)).mean(0)
    if p0:
        for k in range(3): F+=p0[0]*PI[k]*norm.cdf((right-MU[k])/np.sqrt(SD[k]**2+p0[1]))
    return F
def W1pair(a,b): return np.abs(np.cumsum(a)-np.cumsum(b)).sum()*dx
def sample_prehLaw(rng,r,t_prev,h,init,n):
    comps,p0=Kcomp(t_prev,h,init); ws=np.array([w for w,v in comps]+([p0[0]] if p0 else [])); vs=np.array([v for w,v in comps]+([p0[1]] if p0 else []))
    j=rng.choice(len(ws),size=n,p=ws/ws.sum()); base=r[rng.integers(0,len(r),n)]
    if p0: isP0=(j==len(ws)-1); base=np.where(isP0,draw_P0(rng,(n,))[0],base)
    noise=np.where(rng.random(n)<LAM,0.0,np.sqrt(vs[j])*rng.standard_normal(n))
    return base+noise
