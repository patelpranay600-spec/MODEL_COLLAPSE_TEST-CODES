# Is the finite-N discrepancy W1(B,C) explained by an effective sample size N_eff = N*lam*(2-lam) (Wright-Fisher drift with reinjection)?
from kde_common import *
H=0.1; init='ref'
print("N_eff = N*lam*(2-lam) =",N*LAM*(2-LAM),flush=True)
for m in [20,100,400]:
    rx,rk=arch_x[:,:m],arch_k[:,:m]; Rm=[Rmass_of(rx[i]) for i in range(NSEED)]
    rec,_=run2(rx,rk,np.full(NSEED,H),88000+m,init=init,rec_times=(59,)); X,K=rec[59]; rng=np.random.default_rng(99500+m)
    wC=[];wN=[];wE=[]
    for i in range(NSEED):
        Bq=B_grid(Rm[i],T,H,init); wC.append(W1pair(filt(binpts(X[i]),H),Bq))
        for n,lst in [(N,wN),(int(round(N*LAM*(2-LAM))),wE)]:
            Xs=sample_prehLaw(rng,rx[i],T-1,H,init,n); lst.append(W1pair(filt(binpts(Xs),H),Bq))
    f=lambda a:f"{np.mean(a):.4f}±{np.std(a,ddof=1)/np.sqrt(NSEED):.4f}"
    print(f"m={m}: W1(B,C)={f(wC)}   iid floor n=N: {f(wN)}   iid floor n=N_eff=195: {f(wE)}   paired C-floor(N_eff)={np.mean(np.array(wC)-np.array(wE)):+.4f}±{np.std(np.array(wC)-np.array(wE),ddof=1)/np.sqrt(NSEED):.4f}",flush=True)
