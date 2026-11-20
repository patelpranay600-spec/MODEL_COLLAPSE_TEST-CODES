import pickle, numpy as np
res=pickle.load(open('gmm_phase2b_results.pkl','rb'))
IW,IE,IT,IR,IM,IP,IK=range(7)
def ms(a): a=np.asarray(a,float); return a.mean(), a.std(ddof=1)/np.sqrt(len(a))
print("EXPLORATORY N-scaling of pure recursion (R2), m=100, T=60. growth = late-window(51-60) mean minus value at t=1, paired by seed.\n")
print("learner    N        seeds  W1(t=1)   W1(late)   growth W1          |epsR| growth      P(rare) t=1 -> late   slope W1/gen      growth/growth(N=2000)  sqrt(2000/N)")
for L in ['K3','K2','BIC']:
    base=None; rows=[]
    for Nn in [500,2000,20000,200000]:
        a=np.array([r[3] for r in res if r[0]==Nn and r[1]==L])
        if len(a)==0: continue
        w1=a[:,:,IW]; g_=w1[:,50:].mean(1)-w1[:,0]; ge=np.abs(a[:,50:,IE]).mean(1)-np.abs(a[:,0,IE]); t=np.arange(1,61); t=t-t.mean(); sl=(w1*t).sum(1)/(t**2).sum()
        rows.append((Nn,len(a),w1[:,0].mean(),w1[:,50:].mean(),ms(g_),ms(ge),a[:,0,IP].mean(),a[:,50:,IP].mean(),ms(sl)))
    g2000=[r for r in rows if r[0]==2000][0][4][0]
    for r in rows:
        Nn,n,w0,wl,gg,ge,p0,pl,sl=r
        print(f"{L:6s} {Nn:8d}  {n:4d}   {w0:.3f}    {wl:.3f}    {gg[0]:+.3f}±{gg[1]:.3f}    {ge[0]:+.3f}±{ge[1]:.3f}    {p0:.2f} -> {pl:.2f}      {sl[0]:+.5f}±{sl[1]:.5f}   {gg[0]/g2000:6.2f}            {np.sqrt(2000/Nn):.2f}")
    print()
