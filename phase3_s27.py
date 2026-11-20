# Section 27 / Block 2: reference-budget factorial on arm C. Seeds 61-90. h=0.6 everywhere. Wide grid [-80,90], dx=0.05 (T5(lam=.01)=299 has sd~10).
import os,sys
for _v in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS']: os.environ[_v]='1'
import numpy as np, pickle, time
from multiprocessing import Pool
from scipy.stats import norm
from scipy.ndimage import gaussian_filter1d
import kde_common as kc, phase3_arms as pa, phase3_pilot as pp
H=0.6; N=2000; X0_,X1_,DX=-80.,90.,0.05; nG=int((X1_-X0_)/DX); edges=X0_+DX*np.arange(nG+1); right=edges[1:]
F0=sum(kc.PI[k]*norm.cdf(right,kc.MU[k],kc.SD[k]) for k in range(3))
def filt(m,h): return gaussian_filter1d(m,h/DX,mode='constant',truncate=8.0)
def binm(X,w=None): return np.bincount(np.clip(((X-X0_)/DX).astype(int),0,nG-1),weights=w,minlength=nG)/(len(X) if w is None else 1)
def mets(X,h):
    mR=norm.sf((kc.B23-X)/h).mean(); tau=(norm.cdf((-5.5-X)/h)+norm.sf((9.5-X)/h)).mean(); cum=np.cumsum(filt(binm(X),h))
    return dict(W1=np.abs(cum-F0).sum()*DX,rho3=mR,eR=abs(mR-kc.MR0)/kc.MR0,tau=tau,etau=abs(np.log(tau/kc.TAU0)),rv=(X.var()+h*h)/kc.V0),cum
def qstar(r,lam,h):
    g=filt(binm(r),h); q=np.zeros(nG); wj=lam; 
    while wj>1e-7: q+=wj*g; g=filt(g,h); wj*=(1-lam)
    return np.cumsum(q/q.sum())
T5={0.05:59,0.01:299,0.02:149,0.10:29}
CELLS=[('F1',m,0.05) for m in (20,50,100,200,400)]+[('F2',100,l) for l in (0.01,0.02,0.10)]+[('F3',20,0.01),('F3',40,0.02),('F3',200,0.10)]  # (100,.05) shared; (100,.01..) F2 also reused by none
def sub(seed,*tag): return np.random.default_rng(np.random.SeedSequence([seed,*tag]))
def run_cell(seed,r,rl,lam,pihat):
    nr=int(round(lam*N)); ns=N-nr; T=T5[lam]; init=sub(seed,3); kk=init.choice(3,size=N,p=kc.PI); x0=init.normal(kc.MU[kk],kc.SD[kk])[:ns]; l0=kk[:ns]
    rng=sub(seed,2,len(r),int(round(lam*1000))); i=rng.integers(0,len(r),nr); D=np.concatenate([x0,r[i]]); Lb=np.concatenate([l0,rl[i]]); out={}; wins={'R':range(T-4,T+1),'G':range(16,21)}; acc={'R':[],'G':[]}; cumT=None; om5=None
    for t in range(1,T+1):
        j=rng.integers(0,N,ns); Y=D[j]+H*rng.standard_normal(ns); Ly=Lb[j]; i=rng.integers(0,len(r),nr); D=np.concatenate([Y,r[i]]); Lb=np.concatenate([Ly,rl[i]])
        for k,wr in wins.items():
            if t in wr:
                m,cum=mets(D,H); m['om3']=float((Lb==2).mean()); acc[k].append(m)
                if k=='R' and t==T: cumT=cum; om5=m['om3']
    for k in acc: out[k]={q:float(np.mean([d[q] for d in acc[k]])) for q in acc[k][0]}
    out['om3_T5']=om5; out['pihat']=pihat; out['T']=T; out['cumT']=cumT; return out
def job(seed):
    ar=sub(seed,1); k=ar.choice(3,size=pa.MA,p=kc.PI); arch=ar.normal(kc.MU[k],kc.SD[k]); idx=sub(seed,4).choice(pa.MA,400,replace=False); res=dict(seed=seed,cells={})
    for fam,m,lam in CELLS+[('F1',100,0.05)]:
        key=(m,lam)
        if key in res['cells']: continue
        r=arch[idx[:m]]; rl=k[idx[:m]]; pihat=float((rl==2).mean()); c=run_cell(seed,r,rl,lam,pihat)
        rm,cumR=mets(r,H); cq=qstar(r,lam,H); c['Rh']=rm; c['W1_Qref_P0']=float(np.abs(cq-F0).sum()*DX); c['W1_C_Qref']=float(np.abs(np.cumsum(filt(binm(r),H))*0+c.pop('cumT')-cq).sum()*DX); res['cells'][key]=c
    return res
if __name__=="__main__":
    a,b=int(sys.argv[1]),int(sys.argv[2]); fn=sys.argv[3]; t0=time.time(); R=pickle.load(open(fn,'rb')) if os.path.exists(fn) else []; done={r['seed'] for r in R}; todo=[s for s in range(a,b+1) if s not in done]; print('todo',len(todo),flush=True)
    with Pool(int(os.environ.get('NP',2))) as p:
        for r in p.imap_unordered(job,todo,chunksize=1): R.append(r); pickle.dump(R,open(fn,'wb')); print(len(R),round(time.time()-t0),'s',flush=True)
    print('done')
