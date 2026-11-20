import os
for _v in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS']: os.environ[_v]='1'
import numpy as np, warnings, sys, time, pickle
from multiprocessing import Pool
from scipy.stats import norm
from sklearn.mixture import GaussianMixture
warnings.filterwarnings("ignore")
MU=np.array([-3.,2.,8.]); SD=np.array([1.,.8,.6]); PI=np.array([.70,.27,.03]); B23=5.58249
N=2000; T=60; LAM=0.05; NR=int(LAM*N); NS=N-NR; MS=[20,100,400]; NSEED=60; BASE=2026110800
x0,x1,dx=-20.,30.,0.01; xs=np.arange(x0,x1,dx)+dx/2; right=x0+dx*np.arange(1,len(xs)+1)
F0=sum(PI[k]*norm.cdf(right,MU[k],SD[k]) for k in range(3)); f0=sum(PI[k]*norm.pdf(xs,MU[k],SD[k]) for k in range(3))
TAU0=sum(PI[k]*(norm.cdf(-5.5,MU[k],SD[k])+norm.sf(9.5,MU[k],SD[k])) for k in range(3)); MR0=sum(PI[k]*norm.sf(B23,MU[k],SD[k]) for k in range(3))
V0=0.70*(1+9)+0.27*(0.64+4)+0.03*(0.36+64)-(-1.32)**2
LEARN=['K3','K2','BIC']
def fit(x,L,rs):
    x=x.reshape(-1,1)
    if L=='K3': g=GaussianMixture(3,n_init=1,max_iter=100,random_state=rs).fit(x)
    elif L=='K2': g=GaussianMixture(2,n_init=1,max_iter=100,random_state=rs).fit(x)
    else:
        best=None
        for K in range(1,min(6,max(1,len(x)//3))+1):
            gg=GaussianMixture(K,n_init=1,max_iter=100,random_state=rs).fit(x); b=gg.bic(x)
            if best is None or b<best[0]: best=(b,gg)
        g=best[1]
    return g.weights_.copy(),g.means_.ravel().copy(),np.sqrt(g.covariances_.ravel().copy())
def sample(mod,n,rng):
    w,m,s=mod; k=rng.choice(len(w),size=n,p=w/w.sum()); return rng.normal(m[k],s[k])
def metrics(mod):
    w,m,s=mod; w=w/w.sum()
    F=sum(w[j]*norm.cdf(right,m[j],s[j]) for j in range(len(w))); W1=np.abs(F-F0).sum()*dx
    mR=sum(w[j]*norm.sf(B23,m[j],s[j]) for j in range(len(w))); tau=sum(w[j]*(norm.cdf(-5.5,m[j],s[j])+norm.sf(9.5,m[j],s[j])) for j in range(len(w)))
    mean=(w*m).sum(); var=(w*(s**2+m**2)).sum()-mean**2
    d=sum(w[j]*norm.pdf(xs,m[j],s[j]) for j in range(len(w))); mx=np.where((d[1:-1]>d[:-2])&(d[1:-1]>=d[2:]))[0]+1; mn=np.where((d[1:-1]<d[:-2])&(d[1:-1]<=d[2:]))[0]+1
    edges=np.concatenate([[0],mn,[len(d)]]); basins=sum(1 for a,b in zip(edges[:-1],edges[1:]) if d[a:b].sum()*dx>=0.005 and ((mx>=a)&(mx<b)).any())
    rare=float(((m>B23)&(w>=0.01)).any())
    return np.array([W1,(mR-MR0)/MR0,np.log(tau/TAU0) if tau>0 else np.nan,var/V0,basins,rare,len(w)])
NAMES=['W1','epsR','ltau','rv','modes','rare_present','K']
def job(args):
    m,seed=args; rng=np.random.default_rng(BASE+seed*1000+m); rs0=BASE+seed
    k=rng.choice(3,size=m,p=PI); ref=rng.normal(MU[k],SD[k]); out={}
    kk=rng.choice(3,size=N,p=PI); fresh=rng.normal(MU[kk],SD[kk]); boot=ref[rng.integers(0,m,N)]
    for L in LEARN:
        rs=(rs0*7+LEARN.index(L))%1000000000
        mod0=fit(ref,L,rs); out[(L,'R0')]=metrics(mod0)
        out[(L,'R4')]=metrics(fit(fresh,L,rs)); out[(L,'R5')]=metrics(fit(boot,L,rs))
        for reg in ['R2','R3']:
            mod=mod0; tr=[]
            for t in range(1,T+1):
                D=sample(mod,NS if reg=='R3' else N,rng)
                if reg=='R3': D=np.concatenate([D,ref[rng.integers(0,m,NR)]])
                mod=fit(D,L,(rs+t)%4294967295); tr.append(metrics(mod))
            out[(L,reg)]=np.array(tr)
    return (m,seed,out)
if __name__=="__main__":
    t0=time.time(); jobs=[(m,s) for m in MS for s in range(NSEED)]
    res=[]
    with Pool(2) as p:
        for i,r in enumerate(p.imap(job,jobs,chunksize=1)):
            res.append(r)
            if i%10==0: print(i+1,'/',len(jobs),round(time.time()-t0),'s',flush=True)
    pickle.dump(res,open('gmm_phase2_results.pkl','wb')); print("done",time.time()-t0,flush=True)
