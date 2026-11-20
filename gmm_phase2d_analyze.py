import pickle, numpy as np
import gmm_phase2 as g
res=pickle.load(open('gmm_phase2d_results.pkl','rb')); N=g.NAMES; iW,iE,iR,iT=N.index('W1'),N.index('epsR'),N.index('rare_present'),N.index('ltau')
VARS=['V0','V1','V2','V3']
def late(a,i): return a[50:60,i].mean()
def ms(x): x=np.asarray(x,float); return f"{np.nanmean(x):+.3f}±{np.nanstd(x,ddof=1)/np.sqrt(np.sum(~np.isnan(x))):.3f}"
for m in [100,400]:
    R=[o for mm,s,o in res if mm==m]; print(f"\n=== m={m}, {len(R)} seeds ===")
    print("R0 (gen-0 fit) W1: %.3f"%np.mean([o[('R0',)][0] for o in R]))
    for reg in ['R2','R3']:
        print(f"-- {reg}\n  var  W1(t=1)  W1(late)  growth W1        |epsR|late  P(rare) t1->late  ltau(med late)  late-R0 W1")
        for v in VARS:
            t1=np.array([o[(v,reg)][0,iW] for o in R]); L=np.array([late(o[(v,reg)],iW) for o in R]); e=np.array([abs(late(o[(v,reg)],iE)) for o in R])
            p1=np.mean([o[(v,reg)][0,iR] for o in R]); pl=np.mean([late(o[(v,reg)],iR) for o in R]); lt=np.nanmedian([late(o[(v,reg)],iT) for o in R]); r0=np.array([o[('R0',)][0] for o in R])
            print(f"  {v}   {t1.mean():.3f}   {L.mean():.3f}   {ms(L-t1)}   {e.mean():.3f}   {p1:.2f}->{pl:.2f}   {lt:+.2f}   {ms(L-r0)}")
        print("  paired contrasts on late W1 (a-b):")
        for a,b in [('V1','V0'),('V2','V0'),('V3','V0'),('V3','V1'),('V3','V2')]:
            d=np.array([late(o[(a,reg)],iW)-late(o[(b,reg)],iW) for o in R]); se=d.std(ddof=1)/np.sqrt(len(d)); lo,hi=d.mean()-1.96*se,d.mean()+1.96*se
            tag='better(material)' if hi<-0.05 else 'better' if hi<0 else 'worse' if lo>0 else ('≈' if lo>=-0.05 and hi<=0.05 else 'inconclusive')
            print(f"    {a}-{b}: {d.mean():+.3f} [{lo:+.3f},{hi:+.3f}]  {tag}")
