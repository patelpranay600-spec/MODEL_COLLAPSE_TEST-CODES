import pickle, numpy as np
from gmm_phase2 import NAMES, LEARN, MS, T
res=pickle.load(open('gmm_phase2_results.pkl','rb'))
W=slice(50,60)  # generations 51..60
def get(m,L,reg,key=None):
    arr=np.array([r[2][(L,reg)] for r in res if r[0]==m]); return arr
def ms(a):
    a=np.asarray(a,float); a=a[~np.isnan(a)]; return a.mean(), a.std(ddof=1)/np.sqrt(len(a))
def f(a): m,s=ms(a); return f"{m:+.3f}±{s:.3f}"
def fl(a): m,s=ms(a); return f"{m:.3f}±{s:.3f}"
def slope(y):  # per seed OLS slope over generations 1..T
    t=np.arange(1,T+1); t=t-t.mean(); return (y*t).sum(1)/(t**2).sum()
IW,IE,IT,IR,IM,IP,IK=range(7)
print("Primary window: generations 51-60 (mean of per-generation metrics). Mean ± se over 60 paired seeds.\n")
for m in MS:
    print(f"==================== m={m} ====================")
    for L in LEARN:
        r0=get(m,L,'R0'); r4=get(m,L,'R4'); r5=get(m,L,'R5'); r2=get(m,L,'R2'); r3=get(m,L,'R3')
        late=lambda a,i: a[:,W,i].mean(1)
        print(f"-- learner {L}")
        print("   regime        W1            epsR(signed)  |epsR|        ltau(med)  rv            modes   P(rare comp)  K")
        for nm,a,static in [('R0 ref-only',r0,True),('R4 fresh P0',r4,True),('R5 ref-boot',r5,True),('R2 t=1',r2[:,0,:],True),('R2 late',r2,False),('R3 late',r3,False)]:
            g=(lambda i: a[:,i]) if static else (lambda i: late(a,i))
            print(f"   {nm:12s} {fl(g(IW))} {f(g(IE))} {fl(np.abs(g(IE)))} {np.nanmedian(g(IT)):+8.3f}  {fl(g(IR))} {fl(g(IM))} {fl(g(IP))} {fl(g(IK))}")
        d=lambda x,y,i: late(x,i)-y[:,i]
        print("   contrasts (paired): ",end="")
        print(f"(a) R2late-R2(t=1) W1 {f(late(r2,IW)-r2[:,0,IW])}, |epsR| {f(np.abs(late(r2,IE))-np.abs(r2[:,0,IE]))}, P(rare) {f(late(r2,IP)-r2[:,0,IP])}")
        print(f"                       (b) slope/gen R2: W1 {f(slope(r2[:,:,IW]))}, |epsR| {f(slope(np.abs(r2[:,:,IE])))}, P(rare) {f(slope(r2[:,:,IP]))};  R3: W1 {f(slope(r3[:,:,IW]))}, |epsR| {f(slope(np.abs(r3[:,:,IE])))}")
        print(f"                       (c) R2late-R5 W1 {f(late(r2,IW)-r5[:,IW])}, R2late-R0 W1 {f(late(r2,IW)-r0[:,IW])}, R2late-R5 |epsR| {f(np.abs(late(r2,IE))-np.abs(r5[:,IE]))}, P(rare) R2late-R5 {f(late(r2,IP)-r5[:,IP])}")
        print(f"                       (d) R3late-R2late W1 {f(late(r3,IW)-late(r2,IW))}, |epsR| {f(np.abs(late(r3,IE))-np.abs(late(r2,IE)))};  (e) R2late-R4 W1 {f(late(r2,IW)-r4[:,IW])}")
        # exit condition
        a_=late(r2,IW)-r2[:,0,IW]; b_=slope(r2[:,:,IW]); c_=late(r2,IW)-r5[:,IW]
        ci=lambda v:(v.mean()-1.96*v.std(ddof=1)/np.sqrt(len(v)),v.mean()+1.96*v.std(ddof=1)/np.sqrt(len(v)))
        ae=late(r2,IE); aE=np.abs(ae)-np.abs(r2[:,0,IE]); bE=slope(np.abs(r2[:,:,IE])); cE=np.abs(ae)-np.abs(r5[:,IE])
        okW=ci(a_)[0]>0 and ci(b_)[0]>0 and ci(c_)[0]>0; okE=ci(aE)[0]>0 and ci(bE)[0]>0 and ci(cE)[0]>0
        print(f"   EXIT CONDITION: W1 route {'MET' if okW else 'not met'}; rare-region-mass route {'MET' if okE else 'not met'}\n")
    # trajectories (R2 pure recursion W1 and P(rare), seed-mean)
    print("   R2 trajectories (seed means) at t=1,5,10,20,30,40,50,60:")
    for L in LEARN:
        r2=get(m,L,'R2'); ts=[0,4,9,19,29,39,49,59]
        print(f"   {L}: W1 "+" ".join(f"{r2[:,t,IW].mean():.3f}" for t in ts)+" | P(rare) "+" ".join(f"{r2[:,t,IP].mean():.2f}" for t in ts)+" | modes "+" ".join(f"{r2[:,t,IM].mean():.2f}" for t in ts))
    print()
