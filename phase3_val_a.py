# V18-V20: S-NP score validation. seed namespace 90001+
import numpy as np
from scipy.special import logsumexp
from scipy.stats import norm
MU=np.array([-3.,2.,8.]); SD=np.array([1.,.8,.6]); PI=np.array([.70,.27,.03])
def ref(m,seed):
    r=np.random.default_rng(90001+seed*131+m); k=r.choice(3,size=m,p=PI); return r.normal(MU[k],SD[k]),k
def logp(x,r,s): z=(x[:,None]-r[None,:])/s; return logsumexp(-0.5*z*z,axis=1)-np.log(len(r))-np.log(s*np.sqrt(2*np.pi))
def score(x,r,s):
    z=(r[None,:]-x[:,None]); lw=-0.5*(z/s)**2; w=np.exp(lw-logsumexp(lw,axis=1,keepdims=True)); return (w*z).sum(1)/s**2, (w*r[None,:]).sum(1)
def true_score(x,s):  # score of P0*N_s
    v=SD**2+s**2; lw=np.log(PI)+norm.logpdf(x[:,None],MU,np.sqrt(v)); w=np.exp(lw-logsumexp(lw,axis=1,keepdims=True)); return (w*(-(x[:,None]-MU)/v)).sum(1)
X=np.linspace(-15,25,801); out=[]
print("V18(a,b): analytic vs finite-difference / posterior mean")
worst_a=0; worst_b=0
for m in [20,100,400]:
    for sd in range(5):
        r,_=ref(m,sd)
        for s in [0.1,0.3,1,3,6]:
            h=1e-5; fd=(logp(X+h,r,s)-logp(X-h,r,s))/(2*h); sc,pm=score(X,r,s)
            mask=np.abs(sc)>1e-3; ra=np.max(np.abs(sc-fd)[mask]/np.abs(sc)[mask]) if mask.any() else 0; aa=np.max(np.abs(sc-fd)[~mask]) if (~mask).any() else 0
            worst_a=max(worst_a,ra if ra>1e-6 else 0, aa/1e-8 if aa>1e-8 else 0)
            rb=np.max(np.abs(X+s*s*sc-pm)); worst_b=max(worst_b,rb)
            out.append((m,s,ra,aa,rb))
print(" max rel err (where |score|>1e-3): %.2e ; max abs err small-score: %.2e ; max |D_s - posterior mean|: %.2e"%(max(o[2] for o in out),max(o[3] for o in out),max(o[4] for o in out)))
print("V18(c): DSM stationarity [MC]")
r,_=ref(100,0); rng=np.random.default_rng(90001); s=1.0; n=200000; idx=rng.integers(0,100,n); xi=rng.standard_normal(n); xn=r[idx]+s*xi
phi=lambda x: np.sin(x/3)/(1+0.02*x*x)
sc,_=score(xn,r,s); J0=np.mean((sc+xi/s)**2)
Js=[np.mean((sc+d*phi(xn)+xi/s)**2) for d in [-0.05,-0.02,0.02,0.05]]
print("  J(shat)=%.5f ; perturbed J: %s ; all larger: %s"%(J0,np.round(Js,5),all(j>J0 for j in Js)))
print("V18(d): stability at extreme separation")
xs=np.array([-1e3,0.,1e3]); sc,pm=score(xs,np.array([0.,1.,2.]),0.01); print("  finite:",np.isfinite(sc).all(),np.isfinite(pm).all())
print("V19: s_min by LOO likelihood")
grid=np.exp(np.linspace(np.log(0.05),np.log(6),60))
def loo(r,s):
    z=(r[:,None]-r[None,:])/s; K=np.exp(-0.5*z*z); np.fill_diagonal(K,0); return np.log(K.sum(1)/((len(r)-1)*s*np.sqrt(2*np.pi))+1e-300).sum()
for m in [20,100,400]:
    ss=[];interior=[];sil=[]
    for sd in range(20):
        r,_=ref(m,sd); L=[loo(r,s) for s in grid]; i=int(np.argmax(L)); ss.append(grid[i]); interior.append(0<i<len(grid)-1); sil.append(1.06*r.std(ddof=1)*m**-0.2)
    print("  m=%d: s_min median %.3f [min %.3f, max %.3f]; interior %d/20; Silverman median %.3f"%(m,np.median(ss),min(ss),max(ss),sum(interior),np.median(sil)))
print("V20: S-NP score vs population score of P0*N_s (mean |err| by region, 10 refs); rare-mode score sign at x=8 when ref holds no rare pt")
for m in [20,100,400]:
    for s in [0.3,1.0,3.0]:
        E=np.zeros(3); neg=0; nn=0; sg=[]
        for sd in range(10):
            r,k=ref(m,sd); e=np.abs(score(X,r,s)[0]-true_score(X,s))
            E+=np.array([e[(X<-0.105)&(X>-8)].mean(),e[(X>=-0.105)&(X<5.58)].mean(),e[(X>=5.58)&(X<11)].mean()])/10
            if not (k==2).any(): nn+=1; sg.append(score(np.array([8.0]),r,s)[0][0])
        print("  m=%3d s=%.1f  mean|err| R1 %.3f R2 %.3f R3 %.3f ; no-rare refs %d/10, mean score at x=8: %s ; true score at 8: %.3f"%(m,s,*E,nn,("%.3f"%np.mean(sg)) if sg else "n/a",true_score(np.array([8.0]),s)[0]))
