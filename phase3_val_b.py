import os
for _v in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS']: os.environ[_v]='1'
import numpy as np, pickle, time
from multiprocessing import Pool
from scipy.stats import norm, wasserstein_distance
MU=np.array([-3.,2.,8.]); SD=np.array([1.,.8,.6]); PI=np.array([.70,.27,.03]); B12=-0.10488; B23=5.58249; LAM=0.05; H=0.6
xg=np.arange(-20,30,0.01)+0.005; F0=sum(PI[k]*norm.cdf(xg+0.005,MU[k],SD[k]) for k in range(3)); MR0=sum(PI[k]*norm.sf(B23,MU[k],SD[k]) for k in range(3)); TAU0=sum(PI[k]*(norm.cdf(-5.5,MU[k],SD[k])+norm.sf(9.5,MU[k],SD[k])) for k in range(3))
def ref(m,seed):
    r=np.random.default_rng(90001+seed*131+m); k=r.choice(3,size=m,p=PI); return r.normal(MU[k],SD[k]),k
def score(x,r,s):
    z=r[None,:]-x[:,None]; lw=-0.5*(z/s)**2; lw-=lw.max(1,keepdims=True); w=np.exp(lw); w/=w.sum(1,keepdims=True); return (w*z).sum(1)/s**2
grid=np.exp(np.linspace(np.log(0.05),np.log(6),60))
def loo(r,s):
    z=(r[:,None]-r[None,:])/s; K=np.exp(-0.5*z*z); np.fill_diagonal(K,0); return np.log(K.sum(1)/((len(r)-1)*s*np.sqrt(2*np.pi))+1e-300).sum()
def smin(r): return grid[int(np.argmax([loo(r,s) for s in grid]))]
def correct(y,r,s0,sL,c,K,L,rng):
    x=y+s0*rng.standard_normal(len(y)); lv=np.exp(np.linspace(np.log(s0),np.log(sL),L)) if s0>sL else np.full(L,sL)
    for s in lv:
        a=c*s*s
        for _ in range(K): x=x+0.5*a*score(x,r,s)+np.sqrt(a)*rng.standard_normal(len(x))
    return x
def region(x): return np.where(x<B12,0,np.where(x<B23,1,2))
def w1_P0(x):
    xs=np.sort(x); Fn=np.searchsorted(xs,xg+0.005,side='right')/len(xs); return np.abs(Fn-F0).sum()*0.01
def inputs(kind,n,rng):
    k=rng.choice(3,size=n,p=PI); x=rng.normal(MU[k],SD[k])
    if kind=='Qstar':
        j=np.minimum(rng.geometric(LAM,size=n)-1,400); x=x+np.sqrt((j+1)*H*H)*rng.standard_normal(n)
    return x,k
def job(a):
    m,seed=a; r,_=ref(m,seed); sm=smin(r); rng=np.random.default_rng(90001+500000+seed*97+m); res={'m':m,'seed':seed,'smin':sm,'has_rare':bool((r>B23).any())}
    ref_law=r[rng.integers(0,m,20000)]+sm*rng.standard_normal(20000)
    for kind in ['P0','Qstar']:
        y,k=inputs(kind,2000,rng); ry=region(y); A_in=np.array([[np.mean((k==a)&(ry==b)) for b in range(3)] for a in range(3)])
        res[(kind,'in')]=dict(W1P0=w1_P0(y),mR=np.mean(y>=B23),tau=np.mean((y<=-5.5)|(y>=9.5)),A=A_in)
        for s0n,s0 in [('smin',sm),('1.0',max(1.0,sm)),('2.5',2.5),('6.0',6.0)]:
            for c in [0.05,0.2]:
                o=correct(y,r,s0,sm,c,10,10,rng); ro=region(o); d=np.abs(o-y)
                A=np.array([[np.mean((k==a)&(ro==b)) for b in range(3)] for a in range(3)])
                res[(kind,s0n,c)]=dict(disp=[d[k==a].mean() if (k==a).any() else np.nan for a in range(3)],p95=[np.quantile(d[k==a],.95) if (k==a).any() else np.nan for a in range(3)],A=A,W1P0=w1_P0(o),W1ref=wasserstein_distance(o,ref_law),mR=np.mean(o>=B23),tau=np.mean((o<=-5.5)|(o>=9.5)),var=o.var(),mean=o.mean(),out=o[:0])
    return res
if __name__=="__main__":
    # V21: discretised Langevin stationary variance for a Gaussian score
    print("V21: ULA stationary variance, Gaussian N(0,v), 1e5 chains x 200 steps"); rng=np.random.default_rng(90001+900000)
    for ar in [0.05,0.2,0.5]:
        v=1.7; a=ar*v; x=rng.normal(0,np.sqrt(v),100000)*3
        for _ in range(200): x=x+0.5*a*(-x/v)+np.sqrt(a)*rng.standard_normal(len(x))
        pred=v/(1-a/(4*v)); print("  alpha/v=%.2f  sample var %.4f  predicted %.4f  rel diff %.3f%%"%(ar,x.var(),pred,100*abs(x.var()/pred-1)))
    t0=time.time(); jobs=[(m,s) for m in [20,100,400] for s in range(10)]; R=[]
    with Pool(2) as p:
        for i,r in enumerate(p.imap(job,jobs,chunksize=1)):
            R.append(r); print(i+1,len(jobs),round(time.time()-t0),'s',flush=True)
    pickle.dump(R,open('phase3_val_b_results.pkl','wb')); print('done')
