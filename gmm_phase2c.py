# Phase 2c (EXPLORATORY): is the N-independent part of K3 degradation an EM-optimiser artifact? K3, m=100, N=20000, R2, T=60, 20 seeds, same references as Phase 2.
import os
for _v in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS']: os.environ[_v]='1'
import numpy as np, pickle, time, warnings
from multiprocessing import Pool
from sklearn.mixture import GaussianMixture
import gmm_phase2 as g
warnings.filterwarnings("ignore")
M=100; T=60; NN=20000
def fit_variant(x,var,rs,prev=None):
    x=x.reshape(-1,1)
    if var=='tight': gm=GaussianMixture(3,n_init=3,tol=1e-6,max_iter=500,random_state=rs).fit(x)
    else:  # warm start from previous generation's parameters, tight tolerance
        w,m,s=prev; gm=GaussianMixture(3,tol=1e-6,max_iter=500,weights_init=w/w.sum(),means_init=m.reshape(-1,1),precisions_init=(1/s**2).reshape(-1,1,1),random_state=rs).fit(x)
    return gm.weights_.copy(),gm.means_.ravel().copy(),np.sqrt(gm.covariances_.ravel().copy())
def job(a):
    var,seed=a; rng=np.random.default_rng(g.BASE+seed*1000+M); k=rng.choice(3,size=M,p=g.PI); ref=rng.normal(g.MU[k],g.SD[k])
    rng=np.random.default_rng(g.BASE+700000+seed); rs=(g.BASE+seed)*7%1000000000
    mod=fit_variant(ref,'tight',rs)   # generation-0 fit identical procedure for both variants
    tr=[]
    for t in range(1,T+1):
        mod=fit_variant(g.sample(mod,NN,rng),var,(rs+t)%4294967295,prev=mod); tr.append(g.metrics(mod))
    return (var,seed,np.array(tr))
if __name__=="__main__":
    t0=time.time(); jobs=[(v,s) for v in ['tight','warm'] for s in range(20)]; res=[]
    with Pool(2) as p:
        for i,r in enumerate(p.imap_unordered(job,jobs,chunksize=1)):
            res.append(r)
            if i%5==0: print(i+1,'/',len(jobs),round(time.time()-t0),'s',flush=True)
    pickle.dump(res,open('gmm_phase2c_results.pkl','wb')); print("done",round(time.time()-t0),flush=True)
