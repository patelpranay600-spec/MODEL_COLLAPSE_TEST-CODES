# Section 28 / Block 3: H1 (noise N vs bandwidth H) feasibility/pairing check. Seeds 301-330. T=20, G window 16-20, kappa in {.3,.6,.9}.
import os,sys
for _v in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS']: os.environ[_v]='1'
import numpy as np, pickle, time
from multiprocessing import Pool
import kde_common as kc, phase3_arms as pa
H0=0.6; LAM=0.05; N=2000; NR=100; NS=1900; T=20; WIN=range(16,21); KAPPAS=[0.3,0.6,0.9]; PI3=0.03
def sub(seed,*tag): return np.random.default_rng(np.random.SeedSequence([seed,*tag]))
def setup(seed):
    ar=sub(seed,1); k=ar.choice(3,size=pa.MA,p=kc.PI); arch=ar.normal(kc.MU[k],kc.SD[k]); idx=sub(seed,4).choice(pa.MA,100,replace=False); r=arch[idx]; rl=k[idx]
    init=sub(seed,3); kk=init.choice(3,size=NS,p=kc.PI); x0=init.normal(kc.MU[kk],kc.SD[kk]); i0=init.integers(0,100,NR); return r,rl,x0,kk,i0
def recur(r,rl,x0,kk,i0,h,sig,mode,rng,coupled=None):
    D=np.concatenate([x0,r[i0]]); Lb=np.concatenate([kk,rl[i0]]); ms=[]
    for t in range(1,T+1):
        j=rng.integers(0,N,NS); 
        if mode=='N': e=rng.standard_normal(NS); xi=rng.standard_normal(NS); Y=D[j]+h*e+sig*xi
        elif mode=='H': Y=D[j]+np.sqrt(h*h+sig*sig)*rng.standard_normal(NS)
        elif mode=='Nc': e=rng.standard_normal(NS); xi=rng.standard_normal(NS); Y=D[j]+h*e+sig*xi
        else: e=rng.standard_normal(NS); xi=rng.standard_normal(NS); Y=D[j]+np.sqrt(h*h+sig*sig)*((h*e+sig*xi)/np.sqrt(h*h+sig*sig))   # Hc: coupled epsilon-tilde
        i=rng.integers(0,NR,NR); D=np.concatenate([Y,r[i]]); Lb=np.concatenate([Lb[j],rl[i]])
        if t in WIN:
            m=kc.kde_metrics(D,h); ms.append(dict(W1=m['W1'],ltau=m['ltau'],rv=m['rv'],eR=m['epsR'],etot=((Lb==2).mean()-PI3)/PI3))
    return {k_:float(np.mean([d[k_] for d in ms])) for k_ in ms[0]},D
def job(seed):
    r,rl,x0,kk,i0=setup(seed); res=dict(seed=seed,k={})
    for ik,kap in enumerate(KAPPAS):
        sig=kap*H0/np.sqrt(1-LAM); h=H0*np.sqrt(1-kap**2)
        N_,_=recur(r,rl,x0,kk,i0,h,sig,'N',sub(seed,2,ik,0)); H_,_=recur(r,rl,x0,kk,i0,h,sig,'H',sub(seed,2,ik,1))
        _,Da=recur(r,rl,x0,kk,i0,h,sig,'Nc',sub(seed,5,ik)); _,Db=recur(r,rl,x0,kk,i0,h,sig,'Hc',sub(seed,5,ik))
        res['k'][kap]=dict(N=N_,H=H_,coupled_maxdiff=float(np.abs(Da-Db).max()),h=h,sigma=sig)
    return res
if __name__=="__main__":
    a,b=int(sys.argv[1]),int(sys.argv[2]); fn=sys.argv[3]; t0=time.time(); R=pickle.load(open(fn,'rb')) if os.path.exists(fn) else []; done={r['seed'] for r in R}; todo=[s for s in range(a,b+1) if s not in done]
    with Pool(2) as p:
        for r in p.imap_unordered(job,todo,chunksize=1): R.append(r); pickle.dump(R,open(fn,'wb')); print(len(R),round(time.time()-t0),'s',flush=True)
    print('done')
