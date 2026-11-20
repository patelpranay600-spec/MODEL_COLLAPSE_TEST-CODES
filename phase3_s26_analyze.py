import pickle, numpy as np
from scipy import stats
R=sorted(pickle.load(open('phase3_s26.pkl','rb')),key=lambda r:r['seed']); n=len(R); print('n',n)
g=lambda a,k:np.array([r[a][k] for r in R])
def ci(x,B=5000,seed=0):
    rng=np.random.default_rng(seed); idx=rng.integers(0,len(x),(B,len(x))); b=x[idx].mean(1); return np.percentile(b,[2.5,97.5])
def wilson(k,n,z=1.96):
    p=k/n; return ((p+z*z/(2*n)-z*np.sqrt(p*(1-p)/n+z*z/(4*n*n)))/(1+z*z/n),(p+z*z/(2*n)+z*np.sqrt(p*(1-p)/n+z*z/(4*n*n)))/(1+z*z/n))
print('\nR-window means: W1 rho3 tau rv')
for a in ['C','S','DR','V','DR*','V*','R_h','R_hS']: print(f'{a:5s}',*[f'{g(a,k).mean():.4f}' for k in ['W1','rho3','tau','rv']])
mt={'C':'R_h','DR':'R_h','V':'R_h','DR*':'R_h','V*':'R_h','S':'R_hS'}
pr=[('DR','C'),('V','C'),('S','C'),('DR','R_h'),('V','R_h'),('S','R_hS')]; res=[]
for x,y in pr:
    d=g(x,'W1')-g(y,'W1'); t,p=stats.ttest_1samp(d,0); res.append((x,y,d,p))
order=np.argsort([r[3] for r in res]); adj=np.empty(6); run=0
for rank,i in enumerate(order): run=max(run,(6-rank)*res[i][3]); adj[i]=min(1,run)
print('\nPRIMARY W1 contrasts (paired; Holm over 6): diff, 95% boot CI, p, Holm p')
for (x,y,d,p),a in zip(res,adj): print(f' {x}-{y}: {d.mean():+.4f} {np.round(ci(d),4)} p={p:.4f} holm={a:.4f}')
print('\nTOST vs matched R (margin 0.040, 90% CI): mean, lo, hi, equiv, paired SD')
for x in ['DR','V','S']:
    d=g(x,'W1')-g(mt[x],'W1'); se=d.std(ddof=1)/np.sqrt(n); tq=stats.t.ppf(.95,n-1); lo,hi=d.mean()-tq*se,d.mean()+tq*se; print(f' {x}: {d.mean():+.4f} {lo:+.4f} {hi:+.4f} equiv={lo>-.04 and hi<.04} sd={d.std(ddof=1):.4f}')
print('\nSECONDARY (descriptive): R_hS-R_h, S-R_h')
for nm,d in [('R_hS-R_h',g('R_hS','W1')-g('R_h','W1')),('S-R_h',g('S','W1')-g('R_h','W1'))]: print(f' {nm}: {d.mean():+.4f} {np.round(ci(d),4)}')
print('\nGUARDS (point-estimate: fail if >10% of seeds, i.e. >=4/30; Wilson 95% reported)')
cv=abs(g('C','rv')-1)
for x in ['C','DR','V','S','DR*','V*']:
    m=mt[x]; f={'variance':abs(g(x,'rv')-1)>cv+0.05,'region3':abs(g(x,'rho3')-g(m,'rho3'))>0.015,'tail':g(x,'tau')/g(m,'tau')>1.35}
    print(f' {x:4s} matched={m}: '+'; '.join(f'{k} {v.sum()}/{n} fail {np.round(wilson(v.sum(),n),2)} {"FAIL" if v.mean()>0.10 else "ok"}' for k,v in f.items()))
print('\nDIAGNOSTICS')
for a in ['DR','DR*']: print(f' {a}: mean ESS {np.mean([r[a]["diag"]["ess"] for r in R]):.0f}/9500, clipped frac {np.mean([r[a]["diag"]["clip"] or 0 for r in R]):.4f}')
for a in ['V','V*']: print(f' {a}: mean acceptance {np.mean([r[a]["diag"]["acc"] for r in R]):.4f}, min per-gen {min(r[a]["diag_min_acc"] for r in R):.4f}')
print(' all finite:',all(r[a]['finite'] for r in R for a in ['C','S','DR','V','DR*','V*']))
print('\nANCESTRY A33 (reported only):',{a:round(float(np.nanmean([r[a]['A'][2,2] for r in R if r['has_rare']])),3) for a in ['C','S','DR','V','DR*','V*']})
print('seeds with no rare ref point: uniform',sum(not r['has_rare'] for r in R),'S',sum(not r['has_rare_S'] for r in R))
