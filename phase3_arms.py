# Phase 3 validation step 6: arm implementations (C, S, DR, V, B, R, DR*, V*) for the 1-D KDE replace-style recursion, with invariants and compute counters.
# Seed namespace 93001+. NOT a performance comparison; endpoints are not reported here.
import os
for _v in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS']: os.environ[_v]='1'
import numpy as np, time
from scipy.stats import norm
from sklearn.cluster import KMeans
from sklearn.linear_model import LogisticRegression
from sklearn.model_selection import KFold
MU=np.array([-3.,2.,8.]); SD=np.array([1.,.8,.6]); PI=np.array([.70,.27,.03])
N=2000; LAM=0.05; NR=int(LAM*N); NS=N-NR; H=0.6; MA=10000
def p0pdf(x): return sum(PI[k]*norm.pdf(x,MU[k],SD[k]) for k in range(3))
class Cnt(dict):
    def add(self,k,v): self[k]=self.get(k,0)+v
def kde_pdf(x,D,h): 
    out=np.empty(len(x))
    for i in range(0,len(x),2000): z=(x[i:i+2000,None]-D[None,:])/h; out[i:i+2000]=np.exp(-0.5*z*z).sum(1)/(len(D)*h*np.sqrt(2*np.pi))
    return out
def score(x,r,s):
    z=r[None,:]-x[:,None]; lw=-0.5*(z/s)**2; lw-=lw.max(1,keepdims=True); w=np.exp(lw); w/=w.sum(1,keepdims=True); return (w*z).sum(1)/s**2
grid=np.exp(np.linspace(np.log(0.05),np.log(6),60))
def loo(r,s):
    z=(r[:,None]-r[None,:])/s; K=np.exp(-0.5*z*z); np.fill_diagonal(K,0); return np.log(K.sum(1)/((len(r)-1)*s*np.sqrt(2*np.pi))+1e-300).sum()
def smin_of(r): return grid[int(np.argmax([loo(r,s) for s in grid]))]
def correct(y,r,s0,sL,c,K,L,rng,cnt):
    x=y+s0*rng.standard_normal(len(y)); lv=np.exp(np.linspace(np.log(s0),np.log(sL),L)) if s0>sL else np.full(L,sL)
    for s in lv:
        a=c*s*s
        for _ in range(K): x=x+0.5*a*score(x,r,s)+np.sqrt(a)*rng.standard_normal(len(x)); cnt.add('score_evals',len(x)*len(r)); cnt.add('langevin_steps',len(x))
    return x
def make_archive(seed):
    rng=np.random.default_rng(93001+seed); k=rng.choice(3,size=MA,p=PI); return rng.normal(MU[k],SD[k]),k
def retain_uniform(arch,ak,m,rng): i=rng.choice(MA,m,replace=False); return arch[i],ak[i],np.full(m,1/m)
def retain_strat(arch,ak,m,gamma,rng):
    km=KMeans(3,n_init=3,random_state=int(rng.integers(1e9))).fit(arch.reshape(-1,1)); st=km.labels_; pk=np.bincount(st,minlength=3)/MA
    wk=pk**gamma; mk=np.maximum(1,np.round(m*wk/wk.sum()).astype(int)); 
    while mk.sum()>m: mk[np.argmax(mk)]-=1
    while mk.sum()<m: mk[np.argmax(wk/mk)]+=1
    idx=np.concatenate([rng.choice(np.where(st==j)[0],mk[j],replace=False) for j in range(3)]); w=np.concatenate([np.full(mk[j],pk[j]/mk[j]) for j in range(3)]); return arch[idx],ak[idx],w/w.sum()
def draw_ref(r,rl,w,rng): i=rng.choice(len(r),NR,replace=True,p=w); return r[i],rl[i],i
def draw_synth(X,L,n,rng): j=rng.integers(0,len(X),n); return X[j]+H*rng.standard_normal(n),L[j]
def pick(cx,cl,wts,n,rng): i=rng.choice(len(cx),n,replace=True,p=wts/wts.sum()); return cx[i],cl[i]
def dr_weights(cx,r,w,Dprev,rng,cnt,rmax=20.):
    syn,_=draw_synth(Dprev[0],Dprev[1],10000,rng); cnt.add('candidate_draws',10000)
    f=lambda x: np.column_stack([(x/5)**p for p in range(1,5)]); Xc=np.vstack([f(r),f(syn)]); y=np.r_[np.ones(len(r)),np.zeros(len(syn))]
    pc=np.zeros(len(cx)); 
    for tr,_ in KFold(5,shuffle=True,random_state=int(rng.integers(1e9))).split(Xc):
        clf=LogisticRegression(max_iter=300).fit(Xc[tr],y[tr]); pc+=clf.predict_proba(f(cx))[:,1]/5; cnt.add('classifier_fits',1)
    rh=np.clip(len(syn)/len(r)*pc/(1-pc+1e-12),0,rmax); return rh
def one_arm(arm,seed,T,m=100,s_ex=None):
    rng=np.random.default_rng(93001+10**6+seed); arch,ak=make_archive(seed); rr=np.random.default_rng(93001+2*10**6+seed)
    if arm=='S': r,rl,w=retain_strat(arch,ak,m,1.0,rr)
    else: r,rl,w=retain_uniform(arch,ak,m,rr)
    sm=smin_of(r); cnt=Cnt(); t0=time.time(); rng_arm=np.random.default_rng(93001+3*10**6+seed)   # same stream for all arms: common random numbers at generation 1
    k=rng_arm.choice(3,size=NS,p=PI); X0=np.concatenate([rng_arm.normal(MU[k],SD[k]),[]]); L0=k
    ri,rli,idx=draw_ref(r,rl,w,rng_arm); D=np.concatenate([X0,ri]); Lb=np.concatenate([L0,rli]); inv=dict(n_ok=True,ref_in_set=True,labels_ok=True,Y1=None)
    if arm=='R': return dict(arm=arm,retained=(r,rl,w),smin=sm,cnt=cnt,inv=inv,time=time.time()-t0,final=(r,rl))
    for t in range(1,T+1):
        Y,Ly=draw_synth(D,Lb,NS,rng_arm); cnt.add('kde_draws',NS)
        if t==1: inv['Y1']=Y.copy()
        if arm in('C','S'): St,Ls=Y,Ly
        elif arm=='B': St=correct(Y,r,max(1.0,sm),sm,0.2,10,10,rng_arm,cnt); Ls=Ly
        else:
            cx,cl=draw_synth(D,Lb,5*NS,rng_arm); cnt.add('kde_draws',5*NS); cnt.add('candidate_draws',5*NS)
            if arm=='DR': wt=dr_weights(cx,r,w,(D,Lb),rng_arm,cnt)+1e-12
            elif arm=='DR*': wt=p0pdf(cx)/np.maximum(kde_pdf(cx,D,H),1e-300); cnt.add('kernel_evals',len(cx)*len(D)); cnt.add('oracle_target_evals',len(cx))
            if arm in('DR','DR*'): St,Ls=pick(cx,cl,wt,NS,rng_arm)
            else:   # V and V*: threshold, refill until NS accepted
                if arm=='V': 
                    z=(r[:,None]-r[None,:])/sm; Kk=np.exp(-0.5*z*z)@w-w; loov=np.log(np.maximum((Kk)/((1-w))/(sm*np.sqrt(2*np.pi)),1e-300)); thr=np.quantile(loov,0.01)   # LOO score on reference points
                    vfun=lambda x: np.log(np.maximum((np.exp(-0.5*((x[:,None]-r[None,:])/sm)**2)@w)/(sm*np.sqrt(2*np.pi)),1e-300)); cnt.add('kernel_evals',len(r)**2)
                else:
                    kk=np.random.default_rng(93001+4*10**6+1).choice(3,60000,p=PI); ys=np.random.default_rng(93001+4*10**6+2).normal(MU[kk],SD[kk]); thr=np.quantile(np.log(p0pdf(ys)),0.01); vfun=lambda x: np.log(np.maximum(p0pdf(x),1e-300))
                acc=[];accl=[];need=NS
                while need>0:
                    bx,bl=draw_synth(D,Lb,2*NS,rng_arm); cnt.add('kde_draws',2*NS); cnt.add('candidate_draws',2*NS); v=vfun(bx); cnt.add('kernel_evals',len(bx)*len(r)) if arm=='V' else cnt.add('oracle_target_evals',len(bx)); keep=v>=thr
                    acc.append(bx[keep][:need]); accl.append(bl[keep][:need]); need-=len(acc[-1])
                St,Ls=np.concatenate(acc),np.concatenate(accl)
        ri,rli,idx=draw_ref(r,rl,w,rng_arm); D=np.concatenate([St,ri]); Lb=np.concatenate([Ls,rli])
        inv['n_ok']&=len(D)==N and len(Lb)==N; inv['ref_in_set']&=bool(np.all(np.isin(ri,r))); inv['labels_ok']&=bool(set(np.unique(Lb))<= {0,1,2}); inv['finite']=bool(np.isfinite(D).all())
    return dict(arm=arm,retained=(r,rl,w),smin=sm,cnt=cnt,inv=inv,time=time.time()-t0,final=(D,Lb))
if __name__=="__main__":
    ARMS=['C','S','DR','V','B','DR*','V*','R']; T=5; res={a:[one_arm(a,s,T) for s in range(5)] for a in ARMS}
    print("INVARIANTS (5 seeds, T=5, m=100, N=2000, n_ref=100, synthetic draws per generation = 1900)")
    for a in ARMS:
        r=res[a]; print(f" {a:3s} |D|=N:{all(x['inv']['n_ok'] for x in r)} ref draws from retained set:{all(x['inv']['ref_in_set'] for x in r)} labels in {{0,1,2}}:{all(x['inv']['labels_ok'] for x in r)} finite:{all(x['inv'].get('finite',True) for x in r)}  retained m={len(r[0]['retained'][0])} sum w={r[0]['retained'][2].sum():.6f}")
    same=lambda a,b: all(np.array_equal(x['retained'][0],y['retained'][0]) for x,y in zip(res[a],res[b]))
    print("same retained set as C (uniform): "+", ".join(f"{a}:{same('C',a)}" for a in ['DR','V','B','DR*','V*','R'])+f" ; S differs by design: {not same('C','S')}")
    print("same generation-1 synthetic input Y across C,B,DR,V,S (common random numbers):",{a:bool(np.array_equal(res['C'][0]['inv']['Y1'],res[a][0]['inv']['Y1'])) for a in ['B','S','DR','V','DR*','V*']})
    print("\nCOMPUTE per run (T=5, mean over seeds): wall s | score evals | kernel evals | langevin steps | kde draws | candidate draws | classifier fits | oracle target evals")
    for a in ARMS:
        r=res[a]; g=lambda k: np.mean([x['cnt'].get(k,0) for x in r]); print(f" {a:3s} {np.mean([x['time'] for x in r]):7.2f} | {g('score_evals'):.2e} | {g('kernel_evals'):.2e} | {g('langevin_steps'):.2e} | {g('kde_draws'):.0f} | {g('candidate_draws'):.0f} | {g('classifier_fits'):.0f} | {g('oracle_target_evals'):.0f}")
    print("\ns_min used by V/B (median over seeds): %.3f"%np.median([x['smin'] for x in res['B']]))
