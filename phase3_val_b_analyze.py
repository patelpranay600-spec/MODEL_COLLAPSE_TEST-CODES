import pickle, numpy as np
R=pickle.load(open('phase3_val_b_results.pkl','rb')); S0=['smin','1.0','2.5','6.0']
def mse(v): v=np.asarray(v,float); v=v[~np.isnan(v)]; return v.mean(), v.std(ddof=1)/np.sqrt(len(v))
print("s_min by m:",{m:"%.2f (%.2f-%.2f)"%(np.median([r['smin'] for r in R if r['m']==m]),min(r['smin'] for r in R if r['m']==m),max(r['smin'] for r in R if r['m']==m)) for m in [20,100,400]})
print("refs with no rare point:",{m:sum(not r['has_rare'] for r in R if r['m']==m) for m in [20,100,400]})
MR0=0.03*0+None if False else None
from phase3_val_b import MR0,TAU0
print("P0 region-3 mass %.5f, tail mass %.5f"%(MR0,TAU0))
for m in [20,100,400]:
    RR=[r for r in R if r['m']==m]; print(f"\n===== m={m} ({len(RR)} refs) =====")
    for kind in ['P0','Qstar']:
        i=[r[(kind,'in')] for r in RR]; print(f"-- input {kind}: input W1->P0 {mse([a['W1P0'] for a in i])[0]:.3f}, input region3 mass {np.mean([a['mR'] for a in i]):.4f}")
        print("   s0     c    | disp(L1,L2,L3) mean | W1->P0  W1->ref | R3 mass  tail  | offdiag dA (sum |dA| off-diag, max) | var")
        for s0 in S0:
            for c in [0.05,0.2]:
                x=[r[(kind,s0,c)] for r in RR]; d=np.nanmean([a['disp'] for a in x],0); W=mse([a['W1P0'] for a in x]); Wr=mse([a['W1ref'] for a in x])
                dA=np.mean([a['A']-r[(kind,'in')]['A'] for a,r in zip(x,RR)],0); off=dA-np.diag(np.diag(dA))
                print(f"   {s0:5s} {c:4.2f} | {d[0]:.2f} {d[1]:.2f} {d[2]:.2f} | {W[0]:.3f}±{W[1]:.3f} {Wr[0]:.3f} | {np.mean([a['mR'] for a in x]):.4f} {np.mean([a['tau'] for a in x]):.4f} | {np.abs(off).sum():.4f} {np.abs(off).max():.4f} | {np.mean([a['var'] for a in x]):.1f}")
# conjecture checks
print("\n=== C1: mean displacement (all labels pooled over refs) nondecreasing in s0 (tolerance -0.01) ===")
viol=0; tot=0
for m in [20,100,400]:
    for kind in ['P0','Qstar']:
        for c in [0.05,0.2]:
            RR=[r for r in R if r['m']==m]; v=[np.nanmean([np.nanmean(r[(kind,s0,c)]['disp']) for r in RR]) for s0 in S0]
            for a,b in zip(v[:-1],v[1:]):
                tot+=1; viol+= (b-a)<-0.01
            print(f"  m={m} {kind} c={c}: "+" ".join("%.2f"%x for x in v))
print("  violations %d of %d adjacent pairs"%(viol,tot))
print("\n=== C2: |W1(out_P0) - W1(out_Qstar)| to P0 (criterion <= 0.05), and R3-mass gap ===")
for m in [20,100,400]:
    RR=[r for r in R if r['m']==m]
    for s0 in S0:
        for c in [0.05,0.2]:
            a=np.mean([r[('P0',s0,c)]['W1P0'] for r in RR]); b=np.mean([r[('Qstar',s0,c)]['W1P0'] for r in RR]); ma=np.mean([r[('P0',s0,c)]['mR'] for r in RR]); mb=np.mean([r[('Qstar',s0,c)]['mR'] for r in RR])
            print(f"  m={m} s0={s0:5s} c={c}: |dW1|={abs(a-b):.3f} ({'met' if abs(a-b)<=0.05 else 'not met'})  R3 mass P0-in {ma:.4f} Q*-in {mb:.4f}")
print("\n=== C3: s0 in {smin, 1.0}: off-diagonal |dA| max <= 0.005 ===")
for m in [20,100,400]:
    RR=[r for r in R if r['m']==m]
    for kind in ['P0','Qstar']:
        for s0 in ['smin','1.0']:
            for c in [0.05,0.2]:
                dA=np.mean([r[(kind,s0,c)]['A']-r[(kind,'in')]['A'] for r in RR],0); off=np.abs(dA-np.diag(np.diag(dA))).max()
                print(f"  m={m} {kind} s0={s0} c={c}: max off-diag |dA|={off:.4f} ({'met' if off<=0.005 else 'not met'})")
