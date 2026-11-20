import os
for _v in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS']: os.environ[_v]='1'
import numpy as np, pickle, time, warnings
from multiprocessing import Pool
from sklearn.mixture import GaussianMixture
import gmm_phase2 as g
warnings.filterwarnings("ignore")
BASE2=2026110900; N=2000; T=60; NR=int(0.05*N); NS=N-NR; MS=[100,400]; NSEED=60
VARS=['V0','V1','V2','V3']
def _p(gm): return gm.weights_.copy(),gm.means_.ravel().copy(),np.sqrt(gm.covariances_.ravel().copy())
def cold(x,rs,n): return GaussianMixture(3,n_init=n,tol=1e-6,max_iter=500,random_state=rs).fit(x)
def warm(x,rs,prev):
    w,m,s=prev; return GaussianMixture(3,tol=1e-6,max_iter=500,weights_init=w/w.sum(),means_init=m.reshape(-1,1),precisions_init=(1/s**2).reshape(-1,1,1),random_state=rs).fit(x)
def fitv(x,v,rs,prev):
    x=x.reshape(-1,1)
    if v=='V0': return _p(GaussianMixture(3,n_init=1,max_iter=100,random_state=rs).fit(x))
    if v=='V1': return _p(cold(x,rs,5))
    if v=='V2': return _p(warm(x,rs,prev))
    a=warm(x,rs,prev); b=cold(x,rs,5); return _p(a if a.score(x)>=b.score(x) else b)
def job(a):
    m,seed=a; rng=np.random.default_rng(BASE2+seed*1000+m); k=rng.choice(3,size=m,p=g.PI); ref=rng.normal(g.MU[k],g.SD[k])
    rs=(BASE2+seed)*7%1000000000; mod0=fitv(ref,'V0',rs,None); out={('R0',):g.metrics(mod0)}
    for v in VARS:
        for reg in ['R2','R3']:
            r=np.random.default_rng(BASE2+500000+seed*100+m+VARS.index(v)*7+(reg=='R3')*3); mod=mod0; tr=[]
            for t in range(1,T+1):
                D=g.sample(mod,NS if reg=='R3' else N,r)
                if reg=='R3': D=np.concatenate([D,ref[r.integers(0,m,NR)]])
                mod=fitv(D,v,(rs+t)%4294967295,mod); tr.append(g.metrics(mod))
            out[(v,reg)]=np.array(tr)
    return (m,seed,out)
if __name__=="__main__":
    t0=time.time(); jobs=[(m,s) for m in MS for s in range(NSEED)]; res=[]
    with Pool(2) as p:
        for i,r in enumerate(p.imap(job,jobs,chunksize=1)):
            res.append(r)
            if i%10==0: print(i+1,'/',len(jobs),round(time.time()-t0),'s',flush=True)
    pickle.dump(res,open('gmm_phase2d_results.pkl','wb')); print("done",round(time.time()-t0),flush=True)
