import pickle, numpy as np
from scipy import stats
R=pickle.load(open('phase3_s25_pilot.pkl','rb')); R=sorted(R,key=lambda r:r['seed']); n=len(R); print('n seeds',n)
cl=np.array([r['clamped'] for r in R]); un=~cl; print('clamped',cl.sum(),'unclamped',un.sum())
g=lambda arm,k:np.array([r[arm][k] for r in R])
def ci(x,B=5000,rng=np.random.default_rng(0)):
    b=[np.mean(x[rng.integers(0,len(x),len(x))]) for _ in range(B)]; return np.percentile(b,[2.5,97.5])
def tost(d,delta=0.040,a=0.05):
    n=len(d); m=d.mean(); se=d.std(ddof=1)/np.sqrt(n); t=stats.t.ppf(1-a,n-1); lo,hi=m-t*se,m+t*se; return m,lo,hi,(lo>-delta and hi<delta)
print('\nMeans over seeds (R window): W1, tau, rho3, rv')
for a in ['C','R_h','R_loo','R_match','R_eff','B_orig','B_match']: print(f'{a:8s}',*[f'{g(a,k).mean():.4f}' for k in ['W1','tau','rho3','rv']])
tau_h=g('R_h','tau')
print('\nH-S1: share of B_orig tail excess explained by R_match')
ex=g('B_orig','tau')-tau_h; ex_m=g('R_match','tau')-tau_h; print(' ratio of means',ex_m.mean()/ex.mean(),'CI(boot of ratio)',np.percentile([ex_m[i].mean()/ex[i].mean() for i in (lambda r:[r.integers(0,n,n) for _ in range(5000)])(np.random.default_rng(1))],[2.5,97.5]))
d=g('B_orig','W1')-g('R_match','W1'); print(' B_orig-R_match W1 mean %.4f CI'%d.mean(),ci(d),'corr tail',np.corrcoef(g('B_orig','tau'),g('R_match','tau'))[0,1])
print('\nH-S2a tail guard (unclamped): ratio<=1.35')
rat=g('B_match','tau')/tau_h
for nm,m in [('unclamped',un),('clamped',cl),('all',np.ones(n,bool))]: print(f' {nm}: mean ratio {rat[m].mean():.3f} frac<=1.35 {np.mean(rat[m]<=1.35):.2f}; frac |dtau|<=0.0045 {np.mean(abs(g("B_match","tau")-tau_h)[m]<=0.0045):.2f}')
print(' B_orig ratio unclamped',(g('B_orig','tau')/tau_h)[un].mean(),'frac<=1.35',np.mean((g('B_orig','tau')/tau_h)[un]<=1.35))
print('\nH-S2b rho3 guard |rho3(B_match)-rho3(R_h)|<=0.015 (unclamped)'); dr=abs(g('B_match','rho3')-g('R_h','rho3')); print(' frac',np.mean(dr[un]<=0.015),'mean',dr[un].mean())
print('\nH-S2c TOST W1(B_match) vs R_eff, all seeds, margin 0.040: mean,lo90,hi90,equiv'); print(tost(g('B_match','W1')-g('R_eff','W1')))
print(' secondary vs R_h unclamped:',tost((g('B_match','W1')-g('R_h','W1'))[un]))
print('\nSuperiority checks (paired W1, B_match minus ref):')
for a in ['R_h','R_loo','R_eff']: x=g('B_match','W1')-g(a,'W1'); print(f' vs {a}: {x.mean():+.4f} CI',ci(x))
x=g('B_match','W1')-g('C','W1'); print(' vs C (different learner):',f'{x.mean():+.4f}',ci(x))
print('\nTrajectory W1 (t=1,5,10,20): C, B_orig, B_match, R_h')
for t in [1,5,10,20]: print(t,*[f"{np.mean([r[a]['traj'][t] for r in R]):.3f}" for a in ['C','B_orig','B_match']],f"{g('R_h','W1').mean():.3f}")
print('\nAncestry A (signed conditional, mean): A33 C, B_orig, B_match'); 
for a in ['C','B_orig','B_match']: print(a,np.nanmean([r[a]['A'][2,2] for r in R if r['has_rare']]))
