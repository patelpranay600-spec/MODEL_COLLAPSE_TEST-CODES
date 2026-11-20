# Diagnostic (post hoc): is the V18(a) miss finite-difference rounding? 4th-order stencil with larger step, and high-precision (mpmath) reference at sample points.
import numpy as np, mpmath as mp
from phase3_val_a import ref,logp,score
X=np.linspace(-15,25,801); worst={}
for s in [0.1,0.3,1,3,6]:
    w=0
    for m in [20,100,400]:
        for sd in range(5):
            r,_=ref(m,sd); sc=score(X,r,s)[0]
            for h in [1e-3]:
                f=lambda d: logp(X+d,r,s); fd=(-f(2*h)+8*f(h)-8*f(-h)+f(-2*h))/(12*h)
                mask=np.abs(sc)>1e-3; w=max(w,np.max(np.abs(sc-fd)[mask]/np.abs(sc)[mask]))
    worst[s]=w
print("4th-order stencil h=1e-3: max rel err by s:",{k:"%.1e"%v for k,v in worst.items()})
mp.mp.dps=40; r,_=ref(100,0); s=0.1; errs=[]
for x in [-3.2,0.4,2.0,7.9,12.0]:
    num=sum(mp.e**(-(mp.mpf(x)-mp.mpf(float(ri)))**2/(2*mp.mpf(s)**2))*(mp.mpf(float(ri))-mp.mpf(x)) for ri in r); den=sum(mp.e**(-(mp.mpf(x)-mp.mpf(float(ri)))**2/(2*mp.mpf(s)**2)) for ri in r)
    exact=float(num/den/mp.mpf(s)**2); errs.append(abs(score(np.array([x]),r,s)[0][0]-exact)/max(abs(exact),1e-300))
print("vs 40-digit mpmath at s=0.1: max rel err %.1e"%max(errs))
