# R5-2 full KDE lineage chain (positions + labels, KDE sampling, reference insertion) and
# R5-3 transition-matrix controls (negative: non-transporting arm; positive: artificial transporter)
import numpy as np, sys, time
from scipy.stats import norm
from scipy.ndimage import gaussian_filter1d
rng=np.random.default_rng(20261108+55)
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
N=2000
print("=== R5-2b  larger-sample rerun of the variance check (20000 chains) ===",flush=True)
for lam,T,rec in [(0.10,100,(40,60,80,100)),(0.20,50,(20,30,40,50))]:
    m=int(round(lam*N)); m3=int(round(0.03*m)); p=m3/m
    ref_k=np.array([2]*m3+[0]*(m-m3)); ref_x=np.where(ref_k==2,rng.normal(8,.6,m),rng.normal(-3,1,m))
    a=1-lam; Vf=a*p*(1-p)/(N*lam*(2-lam)+a); Vw=p*(1-p)/(N*lam*(2-lam)+a)
    for rule in ['fixed','wr']:
        chains=20000
        om,_=simulate(chains,T,lam,0.6,N,ref_x,ref_k,rule,record=rec)
        arr=np.stack([om[t] for t in rec],1); v=arr.var(0).mean()
        se=np.std([arr[rng.integers(0,chains,chains)].var(0).mean() for _ in range(100)])
        print(f"lam={lam} rule={rule:5s}: var {v:.4e} (se {se/v*100:.1f}%)  ratio to V_fixed {v/Vf:.3f}  to V_wr {v/Vw:.3f}",flush=True)
print("\n=== R5-3b  weak transporter (sensitivity), lam=0.05 h=0.6 T=59 ===",flush=True)
lam,h,T=0.05,0.6,59; m=100; a=1-lam
ref_k=rng.choice(3,size=m,p=pis); ref_x=rng.normal(mus[ref_k],sds[ref_k])
dx=0.01; x=np.arange(-40,60,dx)
comp=[pis[k]*norm.pdf(x,mus[k],sds[k]) for k in range(3)]
reg=[x<b12,(x>=b12)&(x<b23),x>=b23]
spk=[]
for k in range(3):
    d=np.zeros_like(x); idx=np.round((ref_x[ref_k==k]-x[0])/dx).astype(int); np.add.at(d,idx,1.0/(m*dx)); spk.append(d)
def grid_A(f):
    q=[gaussian_filter1d(a*comp[k]+lam*spk[k],h/dx,mode='constant') for k in range(3)]
    Rh=[gaussian_filter1d(spk[k],h/dx,mode='constant') for k in range(3)]; phi=norm.pdf(x,8,.6)
    for t in range(1,T+1):
        q0=q[0]; tr=(1-f)*q0+f*q0.sum()*dx*phi
        q=[gaussian_filter1d(a*(tr if k==0 else q[k]),h/dx,mode='constant')+lam*Rh[k] for k in range(3)]
    return np.array([[q[k][reg[r]].sum()*dx for r in range(3)] for k in range(3)])
A_sc=grid_A(0.0)
np.set_printoptions(precision=5,suppress=True,linewidth=140)
chains=3000
_,Aneg=simulate(chains,T,lam,h,N,ref_x,ref_k,'fixed',f=0.0)
for f in [0.002,0.0005]:
    pred=grid_A(f)-A_sc
    _,Apos=simulate(chains,T,lam,h,N,ref_x,ref_k,'fixed',f=f)
    d=Apos.mean(0)-Aneg.mean(0); se=np.sqrt(Apos.var(0,ddof=1)/chains+Aneg.var(0,ddof=1)/chains)
    print(f"f={f}: predicted Delta A[label1, R3]={pred[0,2]:.5f}  (= {pred[0,2]/0.03*100:.1f}% of the target rare-region mass 0.03); empirical {d[0,2]:.5f} (se {se[0,2]:.5f})  detection z={d[0,2]/se[0,2]:.1f}; Delta A[label3,R3] empirical {d[2,2]:.5f} (se {se[2,2]:.5f}), z={d[2,2]/se[2,2]:.1f}",flush=True)
    print("   full predicted Delta A:\n",pred,"\n   empirical:\n",d,flush=True)
