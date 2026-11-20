# Phase 3 pilot (stage 1/2). Pilot seeds 1-100 (protocol seeds). m=100, lambda=0.05, N=2000, h=0.6, T=T5=59, R window 55-59. Arms: C, B (schedule grid), R.
import os,sys
for _v in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS']: os.environ[_v]='1'
import numpy as np, pickle, time
from multiprocessing import Pool
import kde_common as kc
import phase3_arms as pa
M=100; T=int(os.environ.get('PILOT_T',59)); WIN=range(55,60) if T==59 else range(T-1,T+1); H=pa.H; N=pa.N; NS=pa.NS; NR=pa.NR
GRID=[(c,L,K) for (L,K) in [(10,10),(20,20)] for c in (0.05,0.2)]
def sub(seed,tag): return np.random.default_rng(np.random.SeedSequence([seed,tag]))
def s0_list(sm):
    cand=[('smin',sm),('1.0',max(sm,1.0)),('2.5',max(sm,2.5))]; out=[]; 
    for n,v in cand:
        if not any(abs(v-u)<=1e-3*max(v,u) for _,u in out): out.append((n,v))
    return out
def region(x): return np.where(x<kc.B12,0,np.where(x<kc.B23,1,2)) if hasattr(kc,'B12') else np.where(x<-0.10488,0,np.where(x<kc.B23,1,2))
def gen_metrics(D,Lb):
    m=kc.kde_metrics(D,H); reg=region(D); A=np.full((3,3),np.nan); w=np.array([(Lb==k).mean() for k in range(3)])
    for k in range(3):
        if (Lb==k).any(): A[k]=[np.mean(reg[Lb==k]==r) for r in range(3)]
    return dict(W1=m['W1'],epsR=m['epsR'],rho3=kc.MR0*(1+m['epsR']),tau=m['tau'],ltau=m['ltau'],rv=m['rv'],A=A,w=w)
def summarise(ms):
    out={k:float(np.mean([d[k] for d in ms])) for k in ['W1','epsR','rho3','tau','ltau','rv']}; out['A']=np.nanmean([d['A'] for d in ms],0) if not np.all(np.isnan([d['A'][2,0] for d in ms])) else np.full((3,3),np.nan); out['w']=np.mean([d['w'] for d in ms],0)
    out['finite']=all(np.isfinite([d['W1'] for d in ms])); return out
def run_arm(seed,arm,cfg,r,rl,w0,sm):
    rng=sub(seed,2); init=sub(seed,3)       # init stream shared by all arms -> same D_0; recursion stream shared at generation 1
    k=init.choice(3,size=NS,p=pa.PI); X0=init.normal(pa.MU[k],pa.SD[k]); ri,rli,_=pa.draw_ref(r,rl,w0,init); D=np.concatenate([X0,ri]); Lb=np.concatenate([k,rli]); cnt=pa.Cnt(); t0=time.time(); ms=[]
    for t in range(1,T+1):
        Y,Ly=pa.draw_synth(D,Lb,NS,rng)
        if arm=='C': St,Ls=Y,Ly
        else: c,L,K,s0=cfg; St=pa.correct(Y,r,s0,sm,c,K,L,rng,cnt); Ls=Ly
        ri,rli,_=pa.draw_ref(r,rl,w0,rng); D=np.concatenate([St,ri]); Lb=np.concatenate([Ls,rli])
        if t in WIN: ms.append(gen_metrics(D,Lb))
    o=summarise(ms); o['time']=time.time()-t0; o['score_evals']=cnt.get('score_evals',0); o['updates_per_corr']=(cfg[1]*cfg[2]) if cfg else 0; return o
def job(seed):
    arch_rng=sub(seed,1); k=arch_rng.choice(3,size=pa.MA,p=pa.PI); arch=arch_rng.normal(pa.MU[k],pa.SD[k]); ak=k
    r,rl,w0=pa.retain_uniform(arch,ak,M,sub(seed,4)); sm=pa.smin_of(r); res={'seed':seed,'smin':sm,'s0':s0_list(sm),'nominal_changed':[(n,v) for n,v in [('1.0',1.0),('2.5',2.5)] if sm>v]}
    res['C']=run_arm(seed,'C',None,r,rl,w0,sm)
    # R: reference-only KDE at h=0.6 and at LOO bandwidth
    rm=kc.kde_metrics(r,H); res['R_h']=dict(W1=rm['W1'],epsR=rm['epsR'],rho3=kc.MR0*(1+rm['epsR']),tau=rm['tau'],ltau=rm['ltau'],rv=rm['rv'])
    hl=pa.smin_of(r); rm=kc.kde_metrics(r,hl); res['R_loo']=dict(W1=rm['W1'],epsR=rm['epsR'],rho3=kc.MR0*(1+rm['epsR']),tau=rm['tau'],ltau=rm['ltau'],rv=rm['rv'],h=hl)
    res['has_rare']=bool((r>kc.B23).any()); res['B']={}
    for (c,L,K) in GRID:
        for n,s0 in res['s0']: res['B'][(n,c,L,K)]=run_arm(seed,'B',(c,L,K,s0),r,rl,w0,sm)
    return res
if __name__=="__main__":
    a,b=int(sys.argv[1]),int(sys.argv[2]); tag=sys.argv[3]; t0=time.time(); fn=f'phase3_pilot_{tag}.pkl'
    R=pickle.load(open(fn,'rb')) if os.path.exists(fn) else []; done={r['seed'] for r in R}; todo=[s for s in range(a,b+1) if s not in done]; print('resuming; done',sorted(done),flush=True)
    with Pool(2) as p:
        for i,r in enumerate(p.imap_unordered(job,todo,chunksize=1)):
            R.append(r); pickle.dump(R,open(f'phase3_pilot_{tag}.pkl','wb')); print(len(R),b-a+1,round(time.time()-t0),'s',flush=True)
    print('done')
