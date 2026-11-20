# Phase 2b (EXPLORATORY, added after seeing Phase 2 results): N-scaling of pure recursion R2. Same references as Phase 2 (m=100, same seeds).
import os
for _v in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS']: os.environ[_v]='1'
import numpy as np, pickle, time, warnings
from multiprocessing import Pool
import gmm_phase2 as g
warnings.filterwarnings("ignore")
M=100; T=60
CFG=[(500,['K3','K2','BIC'],20),(2000,['K3','K2','BIC'],20),(20000,['K3','K2','BIC'],20),(200000,['K3','K2'],10)]
def job(a):
    Nn,L,seed=a; rng=np.random.default_rng(g.BASE+seed*1000+M)
    k=rng.choice(3,size=M,p=g.PI); ref=rng.normal(g.MU[k],g.SD[k])           # identical reference to Phase 2
    rng=np.random.default_rng(g.BASE+500000+seed*100+Nn%97); rs=(g.BASE+seed)*7%1000000000
    mod=g.fit(ref,L,rs+g.LEARN.index(L)); tr=[]
    for t in range(1,T+1):
        mod=g.fit(g.sample(mod,Nn,rng),L,(rs+t)%4294967295); tr.append(g.metrics(mod))
    return (Nn,L,seed,np.array(tr))
if __name__=="__main__":
    t0=time.time(); jobs=[(Nn,L,s) for Nn,Ls,ns in CFG for L in Ls for s in range(ns)]
    jobs.sort(key=lambda j:-j[0]); res=[]
    with Pool(2) as p:
        for i,r in enumerate(p.imap_unordered(job,jobs,chunksize=1)):
            res.append(r)
            if i%10==0: print(i+1,'/',len(jobs),round(time.time()-t0),'s',flush=True)
    pickle.dump(res,open('gmm_phase2b_results.pkl','wb')); print("done",round(time.time()-t0),flush=True)
