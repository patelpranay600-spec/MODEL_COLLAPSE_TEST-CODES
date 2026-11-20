# Section 26 / Block 1: baselines C, S, DR, V, DR*, V* vs reference-only R_h, R_h^S. Seeds 31-60. m=100, lam=.05, N=2000, h=.6, T=59, R window 55-59.
import os,sys
for _v in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS']: os.environ[_v]='1'
import numpy as np, pickle, time
from multiprocessing import Pool
import kde_common as kc, phase3_arms as pa, phase3_pilot as pp
H=0.6; T=int(os.environ.get('S26_T',59)); WIN=range(55,60) if T==59 else range(T-1,T+1); NS=pa.NS; ARMS=['C','S','DR','V','DR*','V*']
def wmetrics(r,w,h):   # weighted reference-only KDE metrics (== kc.kde_metrics for uniform w)
    from scipy.stats import norm
    mR=(w*norm.sf((kc.B23-r)/h)).sum(); tau=(w*(norm.cdf((-5.5-r)/h)+norm.sf((9.5-r)/h))).sum()
    cum=np.cumsum(kc.filt(np.bincount(np.clip(((r-kc.x0)/kc.dx).astype(int),0,kc.nG-1),weights=w,minlength=kc.nG),h)); mu=(w*r).sum(); var=(w*(r-mu)**2).sum()
    return dict(W1=kc.W1(cum),rho3=kc.MR0*(1+(mR-kc.MR0)/kc.MR0),tau=tau,rv=(var+h*h)/kc.V0)
def run(seed,arm,r,rl,w0,sm):
    rng=pp.sub(seed,2); init=pp.sub(seed,3); k=init.choice(3,size=NS,p=pa.PI); X0=init.normal(pa.MU[k],pa.SD[k]); ri,rli,_=pa.draw_ref(r,rl,w0,init); D=np.concatenate([X0,ri]); Lb=np.concatenate([k,rli]); cnt=pa.Cnt(); win=[]; diag=dict(ess=[],clip=[],acc=[]); t0=time.time()
    if arm=='V':
        z=(r[:,None]-r[None,:])/sm; Kk=np.exp(-0.5*z*z)@w0-w0; loov=np.log(np.maximum(Kk/(1-w0)/(sm*np.sqrt(2*np.pi)),1e-300)); thr=np.quantile(loov,0.01); vfun=lambda x: np.log(np.maximum((np.exp(-0.5*((x[:,None]-r[None,:])/sm)**2)@w0)/(sm*np.sqrt(2*np.pi)),1e-300))
    if arm=='V*':
        kk=np.random.default_rng(93001+4*10**6+1).choice(3,60000,p=pa.PI); ys=np.random.default_rng(93001+4*10**6+2).normal(pa.MU[kk],pa.SD[kk]); thr=np.quantile(np.log(pa.p0pdf(ys)),0.01); vfun=lambda x: np.log(np.maximum(pa.p0pdf(x),1e-300))
    for t in range(1,T+1):
        if arm in('C','S'): St,Ls=pa.draw_synth(D,Lb,NS,rng)
        elif arm in('DR','DR*'):
            cx,cl=pa.draw_synth(D,Lb,5*NS,rng)
            if arm=='DR': rh=pa.dr_weights(cx,r,w0,(D,Lb),rng,cnt); diag['clip'].append(float(np.mean(rh>=20-1e-9))); wt=rh+1e-12
            else: wt=pa.p0pdf(cx)/np.maximum(pa.kde_pdf(cx,D,H),1e-300)
            diag['ess'].append(float(wt.sum()**2/(wt**2).sum())); St,Ls=pa.pick(cx,cl,wt,NS,rng)
        else:
            acc=[];accl=[];need=NS;drawn=0;nacc=0
            while need>0:
                bx,bl=pa.draw_synth(D,Lb,2*NS,rng); drawn+=len(bx); keep=vfun(bx)>=thr; nacc+=int(keep.sum()); acc.append(bx[keep][:need]); accl.append(bl[keep][:need]); need-=len(acc[-1])
            diag['acc'].append(nacc/drawn); St,Ls=np.concatenate(acc),np.concatenate(accl)
        ri,rli,_=pa.draw_ref(r,rl,w0,rng); D=np.concatenate([St,ri]); Lb=np.concatenate([Ls,rli])
        if t in WIN: win.append(pp.gen_metrics(D,Lb))
    o={k_:float(np.mean([d[k_] for d in win])) for k_ in ['W1','rho3','tau','rv']}; o['A']=np.nanmean([d['A'] for d in win],0) if not np.isnan([d['A'][2,0] for d in win]).all() else np.full((3,3),np.nan)
    o['finite']=bool(np.isfinite(D).all()); o['diag']={k_:(float(np.mean(v)) if v else None) for k_,v in diag.items()}; o['diag_min_acc']=float(np.min(diag['acc'])) if diag['acc'] else None; o['time']=time.time()-t0; return o
def job(seed):
    ar=pp.sub(seed,1); k=ar.choice(3,size=pa.MA,p=pa.PI); arch=ar.normal(pa.MU[k],pa.SD[k]); r,rl,w0=pa.retain_uniform(arch,k,100,pp.sub(seed,4)); sm=pa.smin_of(r)
    rS,rlS,wS=pa.retain_strat(arch,k,100,1.0,pp.sub(seed,5)); smS=pa.smin_of(rS)
    res=dict(seed=seed,smin=sm,has_rare=bool((r>kc.B23).any()),has_rare_S=bool((rS>kc.B23).any()),n_strata_rare=int((rS>kc.B23).sum()))
    res['R_h']=wmetrics(r,w0,H); res['R_hS']=wmetrics(rS,wS,H)
    for a in ARMS:
        res[a]=run(seed,a,rS,rlS,wS,smS) if a=='S' else run(seed,a,r,rl,w0,sm)
    return res
if __name__=="__main__":
    a,b=int(sys.argv[1]),int(sys.argv[2]); fn=sys.argv[3]; t0=time.time(); R=pickle.load(open(fn,'rb')) if os.path.exists(fn) else []; done={r['seed'] for r in R}; todo=[s for s in range(a,b+1) if s not in done]; print('todo',todo,flush=True)
    with Pool(2) as p:
        for r in p.imap_unordered(job,todo,chunksize=1): R.append(r); pickle.dump(R,open(fn,'wb')); print(len(R),round(time.time()-t0),'s',flush=True)
    print('done')
