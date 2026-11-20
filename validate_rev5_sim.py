# R5-2 full KDE lineage chain (positions + labels, KDE sampling, reference insertion) and
# R5-3 transition-matrix controls (negative: non-transporting arm; positive: artificial transporter)
import numpy as np, sys, time
from scipy.stats import norm
from scipy.ndimage import gaussian_filter1d
rng=np.random.default_rng(20261108+5)
mus=np.array([-3.,2.,8.]); sds=np.array([1.,.8,.6]); pis=np.array([.70,.27,.03]); b12,b23=-0.10488,5.58249
def draw_P0(shape):
    k=rng.choice(3,size=shape,p=pis); return rng.normal(mus[k],sds[k]),k
def regprobs(x,h):
    c1=norm.cdf((b12-x)/h); c2=norm.cdf((b23-x)/h); return np.stack([c1,c2-c1,1-c2],-1)
def simulate(chains,T,lam,h,N,ref_x,ref_k,rule,f=0.0,record=(),chunk=250):
    m=len(ref_x); nR=int(round(lam*N)); nS=N-nR
    om_rec={t:[] for t in record}; A_acc=[]
    for c0 in range(0,chains,chunk):
        C=min(chunk,chains-c0)
        sx,sk=draw_P0((C,nS))
        def refpart():
            if rule=='fixed':
                assert nR==m
                return np.tile(ref_x,(C,1)),np.tile(ref_k,(C,1))
            idx=rng.integers(0,m,(C,nR)); return ref_x[idx],ref_k[idx]
        rx,rk=refpart()
        X=np.concatenate([sx,rx],1); K=np.concatenate([sk,rk],1)
        for t in range(1,T+1):
            idx=rng.integers(0,N,(C,nS))
            y=np.take_along_axis(X,idx,1)+h*rng.standard_normal((C,nS)); yk=np.take_along_axis(K,idx,1)
            if f>0:
                mv=(yk==0)&(rng.random((C,nS))<f)
                y=np.where(mv,rng.normal(8.,.6,(C,nS)),y)
            rx,rk=refpart()
            X=np.concatenate([y,rx],1); K=np.concatenate([yk,rk],1)
            if t in om_rec: om_rec[t].append((K==2).mean(1))
        # empirical transition matrix at T
        RP=regprobs(X,h)  # C,N,3
        A=np.zeros((C,3,3))
        for k in range(3): A[:,k,:]=((K==k)[:,:,None]*RP).sum(1)/N
        A_acc.append(A)
    return {t:np.concatenate(v) for t,v in om_rec.items()}, np.concatenate(A_acc)
# ---------- R5-2
print("=== R5-2  full KDE chain: stationary variance of lineage fraction vs V_fixed and V_wr ===",flush=True)
N=2000
for lam,T,rec,chains in [(0.05,140,(60,80,100,120,140),1500),(0.10,100,(40,60,80,100),1500),(0.20,50,(20,30,40,50),1500)]:
    m=int(round(lam*N)); m3=int(round(0.03*m)); p=m3/m
    ref_k=np.array([2]*m3+[0]*(m-m3)); ref_x=np.where(ref_k==2,rng.normal(8,.6,m),rng.normal(-3,1,m))
    a=1-lam; Vf=a*p*(1-p)/(N*lam*(2-lam)+a); Vw=p*(1-p)/(N*lam*(2-lam)+a)
    for rule in ['fixed','wr']:
        t0=time.time(); om,_=simulate(chains,T,lam,0.6,N,ref_x,ref_k,rule,record=rec)
        arr=np.stack([om[t] for t in rec],1)           # chains x times
        v=arr.var(0).mean(); mean=arr.mean()
        bs=[arr[rng.integers(0,chains,chains)].var(0).mean() for _ in range(200)]
        se=np.std(bs)
        print(f"lam={lam} m={m} p={p:.4f} rule={rule:5s}: mean {mean:.5f}  var {v:.4e} (bootstrap se {se/v*100:.1f}%)  ratio to V_fixed {v/Vf:.3f}  to V_wr {v/Vw:.3f}   [V_fixed {Vf:.4e}, V_wr {Vw:.4e}]  ({time.time()-t0:.0f}s)",flush=True)
# ---------- R5-3
print("\n=== R5-3  transition-matrix controls, lam=0.05 h=0.6 T=59 N=2000 m=100, fixed reference set ===",flush=True)
lam,h,T=0.05,0.6,59; m=100; a=1-lam
ref_k=rng.choice(3,size=m,p=pis); ref_x=rng.normal(mus[ref_k],sds[ref_k])
print("retained set label counts:",np.bincount(ref_k,minlength=3),flush=True)
dx=0.01; x=np.arange(-40,60,dx)
comp=[pis[k]*norm.pdf(x,mus[k],sds[k]) for k in range(3)]
reg=[x<b12,(x>=b12)&(x<b23),x>=b23]
spk=[]
for k in range(3):
    d=np.zeros_like(x); idx=np.round((ref_x[ref_k==k]-x[0])/dx).astype(int); np.add.at(d,idx,1.0/(m*dx)); spk.append(d)
def grid_A(f):
    q=[gaussian_filter1d(a*comp[k]+lam*spk[k],h/dx,mode='constant') for k in range(3)]
    Rh=[gaussian_filter1d(spk[k],h/dx,mode='constant') for k in range(3)]
    phi=norm.pdf(x,8,.6)
    for t in range(1,T+1):
        q0=q[0]; tr=(1-f)*q0+f*q0.sum()*dx*phi
        q=[gaussian_filter1d(a*(tr if k==0 else q[k]),h/dx,mode='constant')+lam*Rh[k] for k in range(3)]
    return np.array([[q[k][reg[r]].sum()*dx for r in range(3)] for k in range(3)])
A_sc=grid_A(0.0); f=0.2; A_pos_pred=grid_A(f)
chains=1500
_,Aneg=simulate(chains,T,lam,h,N,ref_x,ref_k,'fixed',f=0.0)
_,Apos=simulate(chains,T,lam,h,N,ref_x,ref_k,'fixed',f=f)
def zt(emp,pred):
    mu=emp.mean(0); se=emp.std(0,ddof=1)/np.sqrt(len(emp)); return mu,se,(mu-pred)/se
np.set_printoptions(precision=5,suppress=True,linewidth=140)
mu,se,z=zt(Aneg,A_sc)
print("NEGATIVE control (non-transporting arm C): empirical mean A - seed-conditional comparator A_SC (Delta A):\n",mu-A_sc,"\n se:\n",se,"\n z:\n",np.round(z,2),f"\n max|z|={np.abs(z).max():.2f}, sum z^2={np.sum(z**2):.1f} (9 cells)",flush=True)
mu_p,se_p,z_p=zt(Apos,A_pos_pred)
print("POSITIVE control (transporter, f=0.2 of label-1 synthetic points moved to N(8,0.6^2) each generation):")
print(" predicted Delta A (transporter - C), from grid:\n",A_pos_pred-A_sc)
print(" empirical Delta A (transporter chains - SC comparator):\n",mu_p-A_sc)
print(" empirical vs predicted transporter A, z:\n",np.round(z_p,2),f"\n max|z|={np.abs(z_p).max():.2f}, sum z^2={np.sum(z_p**2):.1f}")
zd=(mu_p-A_sc)/se_p
print(" detection z of empirical Delta A against 0 (cells):\n",np.round(zd,1),flush=True)
print(" row sums emp (neg):",Aneg.mean(0).sum(1)," col sums emp (neg):",Aneg.mean(0).sum(0))
