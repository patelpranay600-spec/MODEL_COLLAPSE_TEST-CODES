# Phase 3 validation step 2: lambda=0 fixed-horizon control (T_fix=60), 1-D KDE, N=2000, h=0.6. Seed namespace 91001+ (stream A) and 92001+ (stream B).
import numpy as np, sys
from scipy.stats import norm
import kde_common as kc
N=2000; T=60; H=0.6; C=300
def run(seed):
    rng=np.random.default_rng(seed); X,_=kc.draw_P0(rng,(C,N)); rec={0:X.copy()}
    for t in range(1,T+1):
        idx=rng.integers(0,N,(C,N)); X=np.take_along_axis(X,idx,1)+H*rng.standard_normal((C,N))
        if t in (1,2,5,10,20,30,40,50,60): rec[t]=X.copy()
    return rec
def modes(x,h=H):
    d=kc.filt(kc.binpts(x),h)/kc.dx; mx=np.where((d[1:-1]>d[:-2])&(d[1:-1]>=d[2:]))[0]+1; mn=np.where((d[1:-1]<d[:-2])&(d[1:-1]<=d[2:]))[0]+1
    e=np.concatenate([[0],mn,[len(d)]]); return sum(1 for a,b in zip(e[:-1],e[1:]) if d[a:b].sum()*kc.dx>=0.005 and ((mx>=a)&(mx<b)).any())
# NOTE: protocol convention. D_0 = N fresh P0 draws (lambda=0). D_t = N draws from KDE_{t-1}(h) ; hat P_t = KDE(D_t,h)
# Predictions [A]: population variance of P_t (KDE incl. h^2): v_t = v0+(t+1)h^2 (infinite N). Finite N: u_t=E var0(D_t)=(1-1/N)(u_{t-1}+h^2), u_0=(1-1/N)v0, var(P_t)=u_t+h^2.
v0=kc.V0; u=(1-1/N)*v0; pred_fin={0:u+H*H}; pred_pop={0:v0+H*H}
for t in range(1,T+1):
    u=(1-1/N)*(u+H*H); pred_fin[t]=u+H*H; pred_pop[t]=v0+(t+1)*H*H
RA=run(91001); RB=run(92001); out=[]
print("t   | var(P_t) MC (A) [se]  MC (B) | pred finite-N | pred population | MC-pred_fin (A) in se | MC-pred_pop")
for t in (0,1,2,5,10,20,30,40,50,60):
    vA=np.array([x.var()+H*H for x in RA[t]]); vB=np.array([x.var()+H*H for x in RB[t]]); se=vA.std(ddof=1)/np.sqrt(C)
    print(f"{t:3d} | {vA.mean():8.3f} [{se:.3f}]  {vB.mean():8.3f} | {pred_fin[t]:8.3f} | {pred_pop[t]:8.3f} | {(vA.mean()-pred_fin[t])/se:+.2f} | {vA.mean()-pred_pop[t]:+.3f}")
print("\nendpoints along the fixed horizon (mean over %d chains, stream A / stream B; KDE metrics at h=%.1f)"%(C,H))
print("t   | W1 A  W1 B | R3 mass A  B (P0 %.4f)| tail A  B (P0 %.5f)| P(no pt in R3) A  B | modes A  B"%(kc.MR0,kc.TAU0))
for t in (1,5,10,20,40,60):
    m={}
    for nm,R in (('A',RA),('B',RB)):
        ms=[kc.kde_metrics(x,H) for x in R[t]]; m[nm]=(np.mean([a['W1'] for a in ms]),np.mean([a['epsR'] for a in ms]),np.mean([a['tau'] for a in ms]),np.mean([(x>=kc.B23).sum()==0 for x in R[t]]),np.mean([modes(x) for x in R[t][:60]]))
    print(f"{t:3d} | {m['A'][0]:.3f} {m['B'][0]:.3f} | {kc.MR0*(1+m['A'][1]):.4f} {kc.MR0*(1+m['B'][1]):.4f} | {m['A'][2]:.5f} {m['B'][2]:.5f} | {m['A'][3]:.3f} {m['B'][3]:.3f} | {m['A'][4]:.2f} {m['B'][4]:.2f}")
# W1 stream agreement: per-t difference of means in se units
print("\nstream reproducibility (W1 at t=60 per-chain): A %.3f ± %.3f ; B %.3f ± %.3f"%tuple(x for R in (RA,RB) for x in (np.mean([kc.kde_metrics(c,H)['W1'] for c in R[60]]),np.std([kc.kde_metrics(c,H)['W1'] for c in R[60]],ddof=1)/np.sqrt(C))))
