from kde_common import *
print("=== PHASE 1, EXPERIMENT 2: lineage-resolved finite-N analysis. h=0.1, protocol initial law, with-replacement reference rule, T=59, 400 seeds (seed 77003) ===",flush=True)
H=0.1; S=400; ax,ak=draw_P0(np.random.default_rng(77003),(S,ARCH))
P03=np.diff(np.concatenate([[0.03*norm.cdf(edges[0],8,.6)],0.03*norm.cdf(right,8,.6)]))
def rho_pop(R3):
    q=filt(A_*P03+LAM*R3,H); Rh=filt(R3,H)
    for _ in range(T): q=filt(A_*q,H)+LAM*Rh
    cum=np.cumsum(q); return (q.sum()-(np.interp(B23,right,cum)))/q.sum(), q.sum()
for m in [20,100,400]:
    rx,rk=ax[:,:m],ak[:,:m]; rec,om=run2(rx,rk,np.full(S,H),88500+m,init='P0',rec_times=(59,),track=True); X,K=rec[59]
    pihat=np.array([(rk[i]==2).mean() for i in range(S)]); w59=om[:,-1]; everlost=(om[:,1:]==0).any(1)
    cat=np.where(pihat==0,1,np.where(w59==0,2,3))
    rows=[]
    for i in range(S):
        Rall=Rmass_of(rx[i]); R3=np.bincount(np.clip(((rx[i][rk[i]==2]-x0)/dx).astype(int),0,nG-1),minlength=nG)/m
        pm=pop_metrics(A_grid(Rall,T,H,'P0')); rp,mass3=rho_pop(R3)
        cm=kde_metrics(X[i],H); l3=(K[i]==2)
        A33=norm.sf((B23-X[i][l3])/H).sum()/N if l3.any() else 0.0
        rows.append(dict(pm=pm,cm=cm,rp=rp,mass3=mass3,A33=A33,w=w59[i],nl3=l3.sum(),out=(X[i][l3]<B23).mean() if l3.any() else np.nan))
    print(f"\n--- m={m}: seeds={S}; P(no rare point in retained set)={np.mean(cat==1):.3f}; P(rare in set, lineage absent at t=59)={np.mean(cat==2):.3f}; P(rare in set, ever absent at some t)={np.mean((pihat>0)&everlost):.3f}; P(survives)={np.mean(cat==3):.3f}",flush=True)
    print("  category                      n    mean w59 (pop pred)  | W1 pop   W1 C   | epsR pop   epsR C   | ltau pop  ltau C(median) | survivors: own-region retention rho3 emp vs pop (diff±se)   frac label-3 points outside R3",flush=True)
    for c,name in [(1,'(i) no rare lineage in set'),(2,'(ii) rare in set, lost by t=59'),(3,'(iii) rare lineage survives')]:
        idx=[i for i in range(S) if cat[i]==c]
        if not idx: print(f"  {name:30s} 0"); continue
        g=lambda f: np.array([f(rows[i]) for i in idx])
        wp=g(lambda r:r['mass3']); ws=g(lambda r:r['w'])
        s=f"  {name:30s} {len(idx):4d}  {ws.mean():.4f} ({wp.mean():.4f})  | {g(lambda r:r['pm']['W1']).mean():.3f}   {g(lambda r:r['cm']['W1']).mean():.3f}  | {g(lambda r:r['pm']['epsR']).mean():+.3f}    {g(lambda r:r['cm']['epsR']).mean():+.3f}   | {np.nanmean(g(lambda r:r['pm']['ltau'])):+.2f}    {np.nanmedian(g(lambda r:r['cm']['ltau'])):+.2f}"
        if c==3:
            rc=np.array([rows[i]['A33']/rows[i]['w'] for i in idx]); rp_=g(lambda r:r['rp']); dd=rc-rp_
            s+=f"   | {rc.mean():.3f} vs {rp_.mean():.3f}  ({dd.mean():+.4f}±{dd.std(ddof=1)/np.sqrt(len(dd)):.4f})  | {np.nanmean(g(lambda r:r['out'])):.3f}"
        print(s,flush=True)
    # overall excess
    dW=np.array([r['cm']['W1']-r['pm']['W1'] for r in rows]); print(f"  all seeds: finite-N excess W1(C)-W1(pop SC) = {dW.mean():+.4f}±{dW.std(ddof=1)/np.sqrt(S):.4f};  by category: "+", ".join(f"{c}: {dW[cat==c].mean():+.4f}" for c in (1,2,3) if (cat==c).any()),flush=True)
print("DONE",flush=True)
