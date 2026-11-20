import pickle, numpy as np
res=pickle.load(open('gmm_phase2c_results.pkl','rb')); IW,IE,IT,IR,IM,IP,IK=range(7)
ms=lambda a:(np.mean(a),np.std(a,ddof=1)/np.sqrt(len(a)))
b=pickle.load(open('gmm_phase2b_results.pkl','rb'))
print("K3, m=100, N=20000, pure recursion R2, T=60, 20 seeds (same references). growth = late(51-60) - t=1.")
def rep(name,arrs):
    a=np.array(arrs); w=a[:,:,IW]; g=w[:,50:].mean(1)-w[:,0]; ge=np.abs(a[:,50:,IE]).mean(1)-np.abs(a[:,0,IE])
    print(f"{name:34s} W1 growth {ms(g)[0]:+.3f}±{ms(g)[1]:.3f}   |epsR| growth {ms(ge)[0]:+.3f}±{ms(ge)[1]:.3f}   P(rare) t=1 {a[:,0,IP].mean():.2f} -> late {a[:,50:,IP].mean():.2f}   late ltau median {np.median(a[:,50:,IT].mean(1)):+.2f}  late modes {a[:,50:,IM].mean():.2f}")
rep("default EM (Phase 2b, n_init=1)",[r[3] for r in b if r[0]==20000 and r[1]=='K3'])
rep("tight EM (n_init=3, tol=1e-6)",[r[2] for r in res if r[0]=='tight'])
rep("warm-started EM (tol=1e-6)",[r[2] for r in res if r[0]=='warm'])
