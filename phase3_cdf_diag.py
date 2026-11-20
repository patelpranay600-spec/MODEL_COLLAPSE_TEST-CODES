# DIAGNOSTIC (not confirmatory): per-generation CDF / regional W1 for C, R (reference-only), B(s0=max(smin,1.0), c=0.2, 20x20) on pilot seeds 1-30, same streams as phase3_pilot.py.
import os,sys
for _v in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS']: os.environ[_v]='1'
import numpy as np, pickle, time
from multiprocessing import Pool
import kde_common as kc, phase3_arms as pa, phase3_pilot as pp
H=pp.H; T=59; B12=-0.10488; B23=kc.B23; SEL=[1,5,10,20,40,59]
reg=np.where(kc.cent<B12,0,np.where(kc.cent<B23,1,2))
def regional(cum):
    d=cum-kc.F0; ab=np.abs(d)*kc.dx; sg=d*kc.dx; return np.array([ab[reg==k].sum() for k in range(3)]),np.array([sg[reg==k].sum() for k in range(3)])
def law_metrics(cum,X,h):
    r,s=regional(cum); m=kc.kde_metrics(X,h); return dict(reg=r,sreg=s,W1=r.sum(),tau=m['tau'],rv=m['rv'],rho3=kc.MR0*(1+m['epsR']))
def run(seed,arm,cfg,r,rl,w0,sm):
    rng=pp.sub(seed,2); init=pp.sub(seed,3); k=init.choice(3,size=pa.NS,p=pa.PI); X0=init.normal(pa.MU[k],pa.SD[k]); ri,rli,_=pa.draw_ref(r,rl,w0,init); D=np.concatenate([X0,ri]); Lb=np.concatenate([k,rli]); cnt=pa.Cnt(); out=[]; cdfs={}
    for t in range(1,T+1):
        Y,Ly=pa.draw_synth(D,Lb,pa.NS,rng)
        if arm=='C': St,Ls=Y,Ly
        else: c,L,K,s0=cfg; St=pa.correct(Y,r,s0,sm,c,K,L,rng,cnt); Ls=Ly
        ri,rli,_=pa.draw_ref(r,rl,w0,rng); D=np.concatenate([St,ri]); Lb=np.concatenate([Ls,rli])
        cum=np.cumsum(kc.filt(kc.binpts(D),H)); out.append(law_metrics(cum,D,H))
        if t in SEL: cdfs[t]=cum[::5].copy()
    return out,cdfs
def job(seed):
    ar=pp.sub(seed,1); k=ar.choice(3,size=pa.MA,p=pa.PI); arch=ar.normal(pa.MU[k],pa.SD[k]); r,rl,w0=pa.retain_uniform(arch,k,pp.M,pp.sub(seed,4)); sm=pa.smin_of(r); res={'seed':seed,'smin':sm}
    res['C'],res['C_cdf']=run(seed,'C',None,r,rl,w0,sm); res['B'],res['B_cdf']=run(seed,'B',(0.2,20,20,max(sm,1.0)),r,rl,w0,sm)
    for nm,h in [('R_h',H),('R_match',np.sqrt(H*H+sm*sm)),('R_loo',sm)]:
        cum=np.cumsum(kc.filt(kc.binpts(r),h)); res[nm]=law_metrics(cum,r,h); res[nm+'_cdf']=cum[::5].copy()
    return res
if __name__=="__main__":
    t0=time.time(); R=[]
    with Pool(2) as p:
        for i,r in enumerate(p.imap_unordered(job,range(1,31),chunksize=1)):
            R.append(r); pickle.dump(R,open('phase3_cdf_diag.pkl','wb')); print(len(R),30,round(time.time()-t0),'s',flush=True)
    print('done')
