from kde_common import *
print("=== PHASE 1, EXPERIMENT 1: Gaussian scale-mixture equivalence. h=0.1, N=2000, lam=0.05, T=59; 40 paired seeds (archive seed 77001) ===",flush=True)
H=0.1
print("Population identity: q_t = R*K_t with K_t=sum_{j<t} lam a^j N(0,(j+1)h^2) + a^t N(0,(t+1)h^2)   [reference-only initial law, init='ref']")
print("                      q_t = a^{t+1} P0*N(0,(t+1)h^2) + R*sum_{j<=t} lam a^j N(0,(j+1)h^2)   [protocol initial law a*P0+lam*R, init='P0']\n",flush=True)
for m in MS:
    rx,rk=arch_x[:,:m],arch_k[:,:m]; Rm=[Rmass_of(rx[i]) for i in range(NSEED)]
    print(f"--- m={m}",flush=True)
    for init in ['ref','P0']:
        d=[];mx=[]; 
        for i in range(NSEED):
            A=A_grid(Rm[i],T,H,init); B=B_grid(Rm[i],T,H,init); d.append(W1pair(A,B)); mx.append(np.abs(np.cumsum(A)-np.cumsum(B)).max())
        dex=[]; ddisc=[]
        for i in range(5):
            A=A_grid(Rm[i],T,H,init); Fe=B_exact_cdf(rx[i],T,H,init); dex.append(np.abs(np.cumsum(A)-Fe).sum()*dx); Bg=B_grid(Rm[i],T,H,init); ddisc.append(np.abs(np.cumsum(Bg)-Fe).sum()*dx)
        d400=[]
        for i in range(10):
            A=A_grid(Rm[i],400,H,init); B=B_grid(Rm[i],400,H,init); d400.append(W1pair(A,B))
        print(f"  init={init:3s}: W1(A,B) grid-vs-grid T=59: max over seeds {max(d):.2e} (mean {np.mean(d):.2e}); max|dCDF| {max(mx):.2e}; T=400 (10 seeds) max {max(d400):.2e}; vs analytic CDF (5 seeds): W1(A,B_exact) {np.mean(dex):.2e}, W1(B_grid,B_exact) {np.mean(ddisc):.2e}",flush=True)
    # finite-N
    for init in ['ref','P0']:
        rec,_=run2(rx,rk,np.full(NSEED,H),88000+m,init=init,rec_times=(59,))
        X,K=rec[59]; rng=np.random.default_rng(99000+m)
        wC=[];wF=[]
        for i in range(NSEED):
            Bq=B_grid(Rm[i],T,H,init); Cmass=filt(binpts(X[i]),H); wC.append(W1pair(Cmass,Bq))
            Xs=sample_prehLaw(rng,rx[i],T-1,H,init,N); wF.append(W1pair(filt(binpts(Xs),H),Bq))
        wC=np.array(wC); wF=np.array(wF)
        print(f"  init={init:3s}: finite-N  W1(B,C)={wC.mean():.4f}±{wC.std(ddof=1)/np.sqrt(NSEED):.4f}   iid-N floor W1(B, KDE_h of N iid draws from the pre-smoothing law of B)={wF.mean():.4f}±{wF.std(ddof=1)/np.sqrt(NSEED):.4f}   paired excess C-floor={(wC-wF).mean():+.4f}±{(wC-wF).std(ddof=1)/np.sqrt(NSEED):.4f}",flush=True)
print("DONE",flush=True)
