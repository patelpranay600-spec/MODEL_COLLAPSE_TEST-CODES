# Follow-up pilot (Section 25.A). Seeds 201-230. Arms: C, R_h, R_loo, R_match, R_eff, B_orig, B_match.
import os
for _v in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS']: os.environ[_v]='1'
import numpy as np, pickle, time, sys
from multiprocessing import Pool
import kde_common as kc, phase3_arms as pa, phase3_pilot as pp
H=0.6; HMIN=0.15; T=59; WIN=range(55,60); SEL=[1,5,10,20]; NS=pa.NS
def draw(X,L,n,rng,h): j=rng.integers(0,len(X),n); return X[j]+h*rng.standard_normal(n),L[j]
def law(D,Lb,h):
    m=kc.kde_metrics(D,h); reg=pp.region(D); A=np.full((3,3),np.nan); w=np.array([(Lb==k).mean() for k in range(3)])
    for k in range(3):
        if (Lb==k).any(): A[k]=[np.mean(reg[Lb==k]==r) for r in range(3)]
    return dict(W1=m['W1'],rho3=kc.MR0*(1+m['epsR']),tau=m['tau'],rv=m['rv'],A=A,w=w)
def run(seed,arm,h_learn,r,rl,w0,sm):
    rng=pp.sub(seed,2); init=pp.sub(seed,3); k=init.choice(3,size=NS,p=pa.PI); X0=init.normal(pa.MU[k],pa.SD[k]); ri,rli,_=pa.draw_ref(r,rl,w0,init); D=np.concatenate([X0,ri]); Lb=np.concatenate([k,rli]); cnt=pa.Cnt(); traj={}; win=[]
    for t in range(1,T+1):
        Y,Ly=draw(D,Lb,NS,rng,h_learn)
        if arm=='C': St,Ls=Y,Ly
        else: St=pa.correct(Y,r,sm,sm,0.2,20,10,rng,cnt); Ls=Ly
        ri,rli,_=pa.draw_ref(r,rl,w0,rng); D=np.concatenate([St,ri]); Lb=np.concatenate([Ls,rli])
        if t in SEL or t in WIN:
            o=law(D,Lb,h_learn)
            if t in SEL: traj[t]=o['W1']
            if t in WIN: win.append(o)
    out={k_:float(np.mean([d[k_] for d in win])) for k_ in ['W1','rho3','tau','rv']}; out['A']=np.nanmean([d['A'] for d in win],0) if not np.isnan([d['A'][2,0] for d in win]).all() else np.full((3,3),np.nan); out['w']=np.mean([d['w'] for d in win],0); out['traj']=traj; out['h']=h_learn; return out
def job(seed):
    ar=pp.sub(seed,1); k=ar.choice(3,size=pa.MA,p=pa.PI); arch=ar.normal(pa.MU[k],pa.SD[k]); r,rl,w0=pa.retain_uniform(arch,k,100,pp.sub(seed,4)); sm=pa.smin_of(r)
    hfit=float(np.sqrt(max(HMIN**2,H**2-sm**2))); heff=float(np.sqrt(hfit**2+sm**2)); res=dict(seed=seed,smin=sm,hfit=hfit,heff=heff,clamped=bool(sm>np.sqrt(H**2-HMIN**2)),has_rare=bool((r>kc.B23).any()))
    for nm,h in [('R_h',H),('R_loo',sm),('R_match',float(np.sqrt(H**2+sm**2))),('R_eff',heff)]:
        m=kc.kde_metrics(r,h); res[nm]=dict(W1=m['W1'],rho3=kc.MR0*(1+m['epsR']),tau=m['tau'],rv=m['rv'],h=h)
    res['C']=run(seed,'C',H,r,rl,w0,sm); res['B_orig']=run(seed,'B',H,r,rl,w0,sm); res['B_match']=run(seed,'B',hfit,r,rl,w0,sm); return res
if __name__=="__main__":
    t0=time.time(); fn='phase3_s25_pilot.pkl'; R=pickle.load(open(fn,'rb')) if os.path.exists(fn) else []; done={r['seed'] for r in R}; todo=[s for s in range(201,231) if s not in done]; print('todo',len(todo),flush=True)
    with Pool(2) as p:
        for r in p.imap_unordered(job,todo,chunksize=1):
            R.append(r); pickle.dump(R,open(fn,'wb')); print(len(R),30,round(time.time()-t0),'s',flush=True)
    print('done')
