import pickle,numpy as np
from scipy.stats import norm
R=pickle.load(open('phase3_pilot_stage1.pkl','rb')); n=len(R); H=0.6
MU=np.array([-3.,2.,8.]); SD=np.array([1.,.8,.6]); PI=np.array([.70,.27,.03])
def wil(k,n,z=1.96):
    p=k/n; d=1+z*z/n; c=(p+z*z/(2*n))/d; h=z*np.sqrt(p*(1-p)/n+z*z/(4*n*n))/d; return max(0,c-h),min(1,c+h)
def tailmass(s): v=np.sqrt(SD**2+s**2); return sum(PI[k]*(norm.cdf(-5.5,MU[k],v[k])+norm.sf(9.5,MU[k],v[k])) for k in range(3))
cfgs=sorted({k for r in R for k in r['B']},key=lambda k:(k[2],k[1],['smin','1.0','2.5'].index(k[0])))
def mse(x): x=np.asarray(x,float); x=x[~np.isnan(x)]; return x.mean(),1.96*x.std(ddof=1)/np.sqrt(len(x)),len(x)
def AJ(o):
    A=o['A'];w=o['w']
    if np.isnan(A[2]).any(): return (np.nan,)*3
    inflow=sum(w[i]*A[i,2] for i in (0,1) if not np.isnan(A[i,2])); return A[2,2],inflow,inflow-w[2]*(1-A[2,2])
print("=== EVIDENCE TABLE (12 configurations; seeds 1-30; mean ± 95% half-width; n = seeds evaluated) ===")
print("config | n | fail var / tail / R3 [Wilson 95%] | W1 B-C | W1 B-R_h | A33 retention (C: .86) | non-rare inflow to R3 (C .026) | s/run | updates")
tab=[]
for c in cfgs:
    S=[r for r in R if c in r['B']]; k=len(S); b=lambda key:[r['B'][c][key] for r in S]
    fv=sum(abs(r['B'][c]['rv']-1)>abs(r['C']['rv']-1)+0.05 for r in S); ft=sum(r['B'][c]['tau']>0.0248 for r in S); f3=sum(not 0.015<=r['B'][c]['rho3']<=0.060 for r in S)
    dC=mse(np.array(b('W1'))-np.array([r['C']['W1'] for r in S])); dR=mse(np.array(b('W1'))-np.array([r['R_h']['W1'] for r in S])); a33=mse([AJ(r['B'][c])[0] for r in S]); inf=mse([AJ(r['B'][c])[1] for r in S])
    print(f"s0={c[0]:4s} c={c[1]:<4} {c[2]}x{c[3]:<2}| {k} | {fv/k:.2f}[{wil(fv,k)[0]:.2f},{wil(fv,k)[1]:.2f}] / {ft/k:.2f}[{wil(ft,k)[0]:.2f},{wil(ft,k)[1]:.2f}] / {f3/k:.2f}[{wil(f3,k)[0]:.2f},{wil(f3,k)[1]:.2f}] | {dC[0]:+.3f}±{dC[1]:.3f} | {dR[0]:+.3f}±{dR[1]:.3f} | {a33[0]:.2f}±{a33[1]:.2f} (n={a33[2]}) | {inf[0]:.3f}±{inf[1]:.3f} | {np.mean(b('time')):.0f} | {c[2]*c[3]}")
print("\n=== Q1 where is the distribution moved? region masses of the sample D_t (R-window mean), B vs C vs target (0.6998, 0.2702, 0.0300) ===")
def rm(o): 
    A=np.nan_to_num(o['A']); return (o['w'][:,None]*A).sum(0)
print("C  :",np.round(np.mean([rm(r['C']) for r in R],0),4))
for c in [('smin',0.2,10,10),('smin',0.2,20,20),('1.0',0.2,20,20),('2.5',0.05,10,10)]:
    S=[r for r in R if c in r['B']]; print(f"B {c}:",np.round(np.mean([rm(r['B'][c]) for r in S],0),4))
print("(W1 cannot be decomposed by region: per-generation CDFs were not stored. Boundary-CDF errors above are the only region-level evidence available without new compute.)")
print("\n=== Q2 tail inflation ===")
print("Hypothesis (post hoc, from the structure of the operator): B's output is drawn at the final correction noise level s_min, then smoothed again by the KDE (h). Predicted law: P0*N_sqrt(h^2+s_min^2) vs reference-only P0*N_h.")
smv=np.array([r['smin'] for r in R]); pred=np.array([tailmass(np.sqrt(H*H+s*s)) for s in smv]); predR=tailmass(H)
print("analytic tail mass: P0*N_0.6 %.4f | per-seed P0*N_sqrt(h^2+s_min^2) mean %.4f (range %.4f-%.4f)"%(predR,pred.mean(),pred.min(),pred.max()))
print("analytic variance ratio sqrt-model: 1+(h^2+s_min^2)/v0-ish -> (v0+h^2+s_min^2)/v0 mean %.3f ; reference-only (v0+h^2)/v0 %.3f"%(np.mean((8.4407+H*H+smv**2)/8.4407),(8.4407+H*H)/8.4407))
print("config | observed B tail mean | per-seed corr(tail_B, double-smoothing prediction) | B/R_h tail ratio mean (95% HW)")
for c in cfgs:
    S=[(i,r) for i,r in enumerate(R) if c in r['B']]; t=np.array([r['B'][c]['tau'] for _,r in S]); pp=np.array([pred[i] for i,_ in S]); rt=np.array([r['B'][c]['tau']/r['R_h']['tau'] for _,r in S]); m=mse(rt)
    print(f"s0={c[0]:4s} c={c[1]:<4} {c[2]}x{c[3]:<2}| {t.mean():.4f} | {np.corrcoef(t,pp)[0,1]:+.2f} | {m[0]:.2f}±{m[1]:.2f} | obs rv {np.mean([r['B'][c]['rv'] for _,r in S]):.3f}")
print("factor main effects on log(tail_B/tail_Rh), mean over seeds (bootstrap 95% over seeds):")
rng=np.random.default_rng(1)
def lr(c,S): return np.array([np.log(r['B'][c]['tau']/r['R_h']['tau']) for r in S if c in r['B']])
def eff(sel_a,sel_b):
    ds=[]
    for _ in range(2000):
        idx=rng.integers(0,n,n); S=[R[i] for i in idx]; a=np.mean([np.mean(lr(c,S)) for c in cfgs if sel_a(c)]); b=np.mean([np.mean(lr(c,S)) for c in cfgs if sel_b(c)]); ds.append(a-b)
    return np.mean(ds),np.percentile(ds,2.5),np.percentile(ds,97.5)
for nm,fa,fb in [("c=0.2 minus c=0.05",lambda c:c[1]==0.2,lambda c:c[1]==0.05),("20x20 minus 10x10",lambda c:c[2]==20,lambda c:c[2]==10),("s0=1.0 minus s0=smin",lambda c:c[0]=='1.0',lambda c:c[0]=='smin'),("s0=2.5 minus s0=smin",lambda c:c[0]=='2.5',lambda c:c[0]=='smin')]:
    e=eff(fa,fb); print(f"  {nm}: {e[0]:+.2f} [{e[1]:+.2f},{e[2]:+.2f}]")
print("\n=== Q3 what is transported? config-level relation between retention, inflow and W1 gap to reference-only ===")
rows=[]
for c in cfgs:
    S=[r for r in R if c in r['B']]; rows.append((np.nanmean([AJ(r['B'][c])[0] for r in S]),np.nanmean([AJ(r['B'][c])[1] for r in S]),np.nanmean([AJ(r['B'][c])[2] for r in S]),np.mean([r['B'][c]['W1']-r['R_h']['W1'] for r in S]),np.mean([r['B'][c]['tau'] for r in S])))
a=np.array(rows); names=['A33','inflow','J3','W1-R_h','tail']
import itertools
print("across the 12 configurations (Pearson, descriptive only):")
for i,j in [(0,3),(1,3),(2,3),(0,4),(1,4)]: print(f"  corr({names[i]}, {names[j]}) = {np.corrcoef(a[:,i],a[:,j])[0,1]:+.2f}")
print("seed-level, within config s0=1.0 c=0.2 20x20: corr(A33, W1-R_h) and corr(inflow, W1-R_h)")
c=('1.0',0.2,20,20); S=[r for r in R if c in r['B'] and not np.isnan(AJ(r['B'][c])[0])]; y=np.array([r['B'][c]['W1']-r['R_h']['W1'] for r in S]); print("  n=%d  %.2f  %.2f"%(len(S),np.corrcoef([AJ(r['B'][c])[0] for r in S],y)[0,1],np.corrcoef([AJ(r['B'][c])[1] for r in S],y)[0,1]))
print("lineage mass in D_t (w3) B vs C: C %.4f ; B(smin,.2,20x20) %.4f ; B(1.0,.2,20x20) %.4f ; reference lineage-3 fraction proxy via R_h rho3 mean %.4f"%(np.mean([r['C']['w'][2] for r in R]),np.mean([r['B'][('smin',0.2,20,20)]['w'][2] for r in R]),np.mean([r['B'][('1.0',0.2,20,20)]['w'][2] for r in R if ('1.0',0.2,20,20) in r['B']]),np.mean([r['R_h']['rho3'] for r in R])))
