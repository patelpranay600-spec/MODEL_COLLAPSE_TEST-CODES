import pickle, numpy as np
from scipy import stats
R=sorted(pickle.load(open('phase3_s28.pkl','rb')),key=lambda r:r['seed']); n=len(R); print('n',n)
DEL={'ltau':0.112,'etot':0.10,'W1':0.040}
print('deterministic coupled check: max |D_N - D_H| over seeds and kappa =',max(v['coupled_maxdiff'] for r in R for v in r['k'].values()),'(tolerance 1e-10)')
for kap in [0.3,0.6,0.9]:
    print(f'\nkappa={kap} (h={R[0]["k"][kap]["h"]:.4f}, sigma={R[0]["k"][kap]["sigma"]:.4f})')
    print(' endpoint   mean(N-H)   s_d(paired)  sd_unpaired  ratio  n=300 need Delta>=0.190 s_d  n=30 need >=0.601 s_d   Delta_eq  feasible300  feasible30  TOST30(90%CI)  equiv30')
    for e in ['ltau','etot','W1','rv','eR']:
        N_=np.array([r['k'][kap]['N'][e] for r in R]); H_=np.array([r['k'][kap]['H'][e] for r in R]); d=N_-H_; sd=d.std(ddof=1); un=np.sqrt(N_.var(ddof=1)+H_.var(ddof=1)); dl=DEL.get(e,{'rv':0.025,'eR':0.05}.get(e))
        se=sd/np.sqrt(n); tq=stats.t.ppf(.95,n-1); lo,hi=d.mean()-tq*se,d.mean()+tq*se
        print(f' {e:6s} {d.mean():+10.4f} {sd:11.4f} {un:12.4f} {sd/un:6.2f} {0.190*sd:12.4f} {0.601*sd:22.4f} {dl:10.3f} {str(dl>=0.190*sd):>11s} {str(dl>=0.601*sd):>11s}  [{lo:+.4f},{hi:+.4f}] {lo>-dl and hi<dl}')
