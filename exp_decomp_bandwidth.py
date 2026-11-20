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
# ---------------- Experiment 1
print("=== EXPERIMENT 1: population vs finite-N decomposition, fixed h=0.1 (N=2000, lam=0.05, T=59, 40 paired seeds, nested reference sets) ===",flush=True)
H1=0.1; A_traj=pop_traj(P0mass,H1)
print("Level A (P0, deterministic operator): "+"  ".join(f"{s}: W1={A_traj[59 if s=='t59' else 20]['W1']:.4f}" for s in ['t20','t59'])+f"  win55-59 W1={np.mean([A_traj[t]['W1'] for t in range(55,60) if t in A_traj]):.4f}  epsR(t59)={A_traj[59]['epsR']:+.4f}  ltau(t59)={A_traj[59]['ltau']:+.4f}  rv(t59)={A_traj[59]['rv']:.4f}",flush=True)
E1={}
for m in MS:
    rx,rk=arch_x[:,:m],arch_k[:,:m]
    B=[pop_traj(binpts_m,H1) for binpts_m in [np.bincount(np.clip(((rx[i]-x0)/dx).astype(int),0,nG-1),minlength=nG)/m for i in range(NSEED)]]
    rec=run_chains(rx,rk,np.full(NSEED,H1),88000+m); Cm=chain_metrics(rec,np.full(NSEED,H1))
    refkde=[kde_metrics(rx[i],H1) for i in range(NSEED)]; raw=[raw_metrics(rx[i]) for i in range(NSEED)]
    E1[m]=(B,Cm,refkde,raw)
    print(f"\n--- m={m}: raw reference P(zero tail points)={np.mean([r['tau']==0 for r in raw]):.3f}; P(retained set has no rare point)={np.mean([(rk[i]==2).sum()==0 for i in range(NSEED)]):.3f}; mean pi_hat={np.mean([(rk[i]==2).mean() for i in range(NSEED)]):.4f}",flush=True)
    print("  spec metric |    A(P0,pop)   B(Pref,pop)   C(Pref,finiteN) |  RefKDE(h)    raw ref   | B-A (se)          C-B (se)",flush=True)
    for spec in ['t20','t59','win']:
        for key in (KEYS if spec=='win' else ['W1','epsR']):
            a=np.full(NSEED,np.mean([A_traj[t][key] for t in ([20] if spec=='t20' else [59] if spec=='t59' else range(55,60)) if t in A_traj]))
            b=aggpop(B,key,spec); c=agg(Cm,key,spec); rk_=np.array([d[key] for d in refkde]); rw=np.array([d[key] for d in raw])
            print(f"  {spec:3s} {key:5s} | {fm(a)}   {fm(b)}   {fm(c)} | {fm(rk_)}  {fm(rw)} | {fms(b-a)}  {fms(c-b)}",flush=True)
    # zero-tail / lineage
    tau_c=np.array([[d['tau'] for d in Cm[t]] for t in range(55,60) if t in Cm]) if False else np.array([d['tau'] for d in Cm[59]])
    print(f"  finite-N C at t=59: P(model tail < 1e-3*tau0)={np.mean(tau_c<1e-3*TAU0):.3f}; P(D_t has zero tail points)={np.mean([d['emp0'] for d in Cm[59]]):.3f}; P(label-3 absent)={1-np.mean([d['surv'] for d in Cm[59]]):.3f}; mean omega_3={np.mean([d['omega'] for d in Cm[59]]):.4f}",flush=True)
    pos=[d['ltau'] for d in Cm[59]]; print(f"  ltau (C,t59) quantiles 10/50/90%: {np.nanpercentile(pos,[10,50,90]).round(3)}",flush=True)
# ---------------- Experiment 2
print("\n\n=== EXPERIMENT 2: bandwidth rules, recursive C vs reference-only KDE at the same h (paired seeds) ===",flush=True)
HG=np.geomspace(0.03,3.0,60)
def loo_h(x):
    D=(x[:,None]-x[None,:]); best,bh=-1e18,None
    for h in HG:
        K=np.exp(-0.5*(D/h)**2)/(h*np.sqrt(2*np.pi)); np.fill_diagonal(K,0); dens=K.sum(1)/(len(x)-1)
        ll=np.log(np.maximum(dens,1e-300)).sum()
        if ll>best: best,bh=ll,h
    return bh
def silverman(x): s=x.std(ddof=1); iqr=np.subtract(*np.percentile(x,[75,25])); return 0.9*min(s,iqr/1.34)*len(x)**(-0.2)
def infl(x,eps=0.05): return np.sqrt(eps*LAM*x.var(ddof=1))
ORH=[0.05,0.075,0.1,0.15,0.2,0.3,0.45,0.6]
# oracle bandwidths selected on separate selection seeds
orac_ref,orac_C={},{}
for m in MS:
    rx,rk=sel_x[:,:m],sel_k[:,:m]
    orac_ref[m]=min(ORH,key=lambda h: np.mean([kde_metrics(rx[i],h)['W1'] for i in range(NSEL)]))
    sc={}
    for h in ORH:
        rec=run_chains(rx,rk,np.full(NSEL,h),89000+m); sc[h]=np.mean(agg(chain_metrics(rec,np.full(NSEL,h)),'W1','win'))
    orac_C[m]=min(sc,key=sc.get)
print("oracle h (selected on 20 separate seeds): reference-only",orac_ref," recursive C",orac_C,flush=True)
RULES=['fixed0.1','silverman','LOO','inflation','oracle_ref','oracle_C']
for m in MS:
    rx,rk=arch_x[:,:m],arch_k[:,:m]
    hs={'fixed0.1':np.full(NSEED,0.1),'silverman':np.array([silverman(rx[i]) for i in range(NSEED)]),'LOO':np.array([loo_h(rx[i]) for i in range(NSEED)]),
        'inflation':np.array([infl(rx[i]) for i in range(NSEED)]),'oracle_ref':np.full(NSEED,orac_ref[m]),'oracle_C':np.full(NSEED,orac_C[m])}
    raw=[raw_metrics(rx[i]) for i in range(NSEED)]
    print(f"\n--- m={m}  (raw unsmoothed reference: W1={fm([r['W1'] for r in raw])} epsR={fm([r['epsR'] for r in raw])} P(zero tail)={np.mean([r['tau']==0 for r in raw]):.3f})",flush=True)
    print("  rule        mean h | W1: Ref     C       C-Ref (se)       | epsR: Ref    C     | ltau: Ref    C     | rv: Ref    C     | C: P(w=0) mean w3 P(tail<1e-3 tau0)",flush=True)
    for r in RULES:
        h=hs[r]; rec=run_chains(rx,rk,h,88000+m); Cm=chain_metrics(rec,h); refk=[kde_metrics(rx[i],h[i]) for i in range(NSEED)]
        g=lambda key: (np.array([d[key] for d in refk]),agg(Cm,key,'win'))
        w=g('W1'); e=g('epsR'); l=g('ltau'); v=g('rv')
        tau59=np.array([d['tau'] for d in Cm[59]])
        print(f"  {r:10s} {h.mean():6.3f} | {fm(w[0])} {fm(w[1])} {fms(w[1]-w[0])} | {fm(e[0])} {fm(e[1])} | {fm(l[0])} {fm(l[1])} | {fm(v[0])} {fm(v[1])} | {1-np.mean([d['surv'] for d in Cm[59]]):.3f} {np.mean([d['omega'] for d in Cm[59]]):.4f} {np.mean(tau59<1e-3*TAU0):.3f}",flush=True)
print("\nDONE",flush=True)
