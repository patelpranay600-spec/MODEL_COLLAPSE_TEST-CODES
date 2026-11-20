from kde_common import *
print("=== PHASE 1, EXPERIMENT 3: recursion-aware bandwidth neighbourhood (same 40 paired seeds as earlier, archive 77001), full endpoint trade-off ===",flush=True)
def infl(x,eps): return np.sqrt(eps*LAM*x.var(ddof=1))
def silverman(x): s=x.std(ddof=1); iqr=np.subtract(*np.percentile(x,[75,25])); return 0.9*min(s,iqr/1.34)*len(x)**(-0.2)
HG=np.geomspace(0.03,3.0,60)
def loo_h(x):
    D=(x[:,None]-x[None,:]); best,bh=-1e18,None
    for h in HG:
        K=np.exp(-0.5*(D/h)**2)/(h*np.sqrt(2*np.pi)); np.fill_diagonal(K,0); ll=np.log(np.maximum(K.sum(1)/(len(x)-1),1e-300)).sum()
        if ll>best: best,bh=ll,h
    return bh
for m in MS:
    rx,rk=arch_x[:,:m],arch_k[:,:m]; Rm=[Rmass_of(rx[i]) for i in range(NSEED)]
    vh=np.array([rx[i].var(ddof=1) for i in range(NSEED)])
    rules={'fixed0.075':np.full(NSEED,0.075),'fixed0.1':np.full(NSEED,0.1)}
    for e in [0.025,0.05,0.1,0.2]: rules[f'infl eps={e}']=np.array([infl(rx[i],e) for i in range(NSEED)])
    print(f"\n--- m={m}  (mean v_hat={vh.mean():.2f}; P(raw ref has zero tail points)={np.mean([((rx[i]<=-5.5)|(rx[i]>=9.5)).sum()==0 for i in range(NSEED)]):.3f})",flush=True)
    # reference-only panel
    panel={}
    panel['raw']=[raw_metrics(rx[i]) for i in range(NSEED)]
    panel['KDE LOO']=[kde_metrics(rx[i],loo_h(rx[i])) for i in range(NSEED)]
    panel['KDE Silverman']=[kde_metrics(rx[i],silverman(rx[i])) for i in range(NSEED)]
    for nm,h in rules.items(): panel['KDE '+nm]=[kde_metrics(rx[i],h[i]) for i in range(NSEED)]
    for nm in ['fixed0.1','infl eps=0.05']:
        hh=rules[nm]; panel['scale-mix (pop limit) @'+nm]=[pop_metrics(A_grid(Rm[i],400,hh[i],'ref')) for i in range(NSEED)]
    print("  REFERENCE-ONLY estimators:      W1            epsR(signed)   ltau          rv",flush=True)
    refW={}
    for nm,L in panel.items():
        g=lambda k: np.array([d[k] for d in L]); refW[nm]=g('W1')
        print(f"    {nm:34s} {fm(g('W1'))}±{ms(g('W1'))[1]:.3f}  {fm(g('epsR'))}   {np.nanmedian(g('ltau')):9.3f}  {fm(g('rv'))}",flush=True)
    print("  RECURSIVE C by bandwidth rule (window 55-59): mean h | W1  C-Ref(same h) | epsR | ltau median (P(tail<1e-3 tau0)) | rv | mean w3 P(w3=0@59)",flush=True)
    Cres={}
    for nm,h in rules.items():
        rec,om=run2(rx,rk,h,88000+m,init='P0',rec_times=(55,56,57,58,59),track=False)
        per={t:[kde_metrics(rec[t][0][i],h[i]) for i in range(NSEED)] for t in rec}
        g=lambda k: np.mean([[d[k] for d in per[t]] for t in per],0); refk=np.array([kde_metrics(rx[i],h[i])['W1'] for i in range(NSEED)])
        tau59=np.array([d['tau'] for d in per[59]]); K59=rec[59][1]; w3=(K59==2).mean(1)
        Cres[nm]=g('W1')
        print(f"    {nm:16s} {h.mean():.3f} | {g('W1').mean():.4f}  {fms(g('W1')-refk)} | {g('epsR').mean():+.4f} | {np.median(g('ltau')):+.3f} ({np.mean(tau59<1e-3*TAU0):.3f}) | {g('rv').mean():.4f} | {w3.mean():.4f} {np.mean(w3==0):.3f}",flush=True)
    bestC=min(Cres,key=lambda k:Cres[k].mean()); bestR=min(refW,key=lambda k:refW[k].mean()); d=Cres[bestC]-refW[bestR]
    print(f"  best prespecified/rule-based C: {bestC} W1={Cres[bestC].mean():.4f}; best reference-only (any listed, incl. oracle-free rules): {bestR} W1={refW[bestR].mean():.4f}; paired gap C-Ref = {d.mean():+.4f}±{d.std(ddof=1)/np.sqrt(NSEED):.4f}",flush=True)
print("DONE",flush=True)
