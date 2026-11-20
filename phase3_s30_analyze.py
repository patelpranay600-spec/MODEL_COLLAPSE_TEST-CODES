import pickle, numpy as np
from scipy import stats
R=sorted(pickle.load(open('phase3_s30.pkl','rb')),key=lambda r:r['seed']); n=len(R); print('n',n,'seeds',R[0]['seed'],'-',R[-1]['seed'])
cm=max(v['coupled_maxdiff'] for r in R for v in r['k'].values()); print('deterministic coupled check max |D_N-D_H| =',cm,'PASS' if cm<=1e-10 else 'FAIL')
DEL={'ltau':0.112,'etot':0.10,'W1':0.040,'rv':0.025,'eR':0.05}
PRIM={0.3:['ltau'],0.6:['ltau','etot'],0.9:['ltau','etot']}   # etot at kappa=0.3 pre-classified inconclusive by design
out={}
for kap in [0.3,0.6,0.9]:
    print(f'\nkappa={kap}')
    print(' endpoint  mean(N-H)  s_d   feasible(0.190 s_d<=Delta)  Delta   90% CI   TOST equivalent   [role]')
    for e in ['ltau','etot','W1','rv','eR']:
        N_=np.array([r['k'][kap]['N'][e] for r in R]); H_=np.array([r['k'][kap]['H'][e] for r in R]); d=N_-H_; sd=d.std(ddof=1); se=sd/np.sqrt(n); tq=stats.t.ppf(.95,n-1); lo,hi=d.mean()-tq*se,d.mean()+tq*se; dl=DEL[e]
        role='primary' if e in PRIM[kap] else ('primary-pre-classified-inconclusive' if (e=='etot' and kap==0.3) else ('secondary(W1, feasible)' if e=='W1' else 'secondary-inconclusive-by-design'))
        eq=lo>-dl and hi<dl; out[(kap,e)]=eq
        print(f' {e:6s} {d.mean():+9.4f} {sd:7.4f} {str(dl>=0.190*sd):>8s} {dl:8.3f} [{lo:+.4f},{hi:+.4f}] {str(eq):>8s}   {role}')
print('\nIU per kappa over primary endpoints {ltau, etot}:')
for kap in [0.3,0.6,0.9]:
    if kap==0.3: print(f' kappa=0.3: not claimable as primary IU (etot pre-classified inconclusive); ltau equivalent = {out[(kap,"ltau")]}')
    else: print(f' kappa={kap}: IU equivalence = {out[(kap,"ltau")] and out[(kap,"etot")]}')
print('H1 overall (IU over kappa=0.6 and 0.9, both endpoints):',all(out[(k,e)] for k in (0.6,0.9) for e in ('ltau','etot')))
