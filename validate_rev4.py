import numpy as np
from scipy.ndimage import gaussian_filter1d
from scipy.stats import norm
rng=np.random.default_rng(20261108)
print("=== R4-A  lineage chain: exact vs approximate stationary variance ===")
def chain(N,lam,p,mode,chains=40000,T=None):
    a=1-lam; nS=int(round(a*N)); nR=int(round(lam*N))
    T=T or int(np.ceil(np.log(1e-6)/np.log(a)))+50
    w=np.full(chains,0.03)
    m3=int(round(p*nR))
    for t in range(T):
        syn=rng.binomial(nS,w)
        ref = m3 if mode=='fixed' else rng.binomial(nR,p,size=chains)
        w=(syn+ref)/N
    return w
for N,lam,p in [(2000,0.05,0.03),(2000,0.10,0.03),(100,0.20,0.05),(100,0.40,0.05)]:
    a=1-lam
    # fixed reference counts (validated rule n_ref=m=lam*N): p must equal m3/nR
    nR=int(round(lam*N)); pf=round(p*nR)/nR
    w=chain(N,lam,pf,'fixed')
    ex=a*pf*(1-pf)/(N*lam*(2-lam)+a); ap=a*pf*(1-pf)/(N*lam*(2-lam))
    print(f"fixed-ref N={N} lam={lam} p={pf:.4f}: mean {w.mean():.5f} (p) var emp {w.var():.4e} exact {ex:.4e} approx {ap:.4e}")
    w=chain(N,lam,p,'resample')
    exr=p*(1-p)/(N*lam*(2-lam)+a)
    print(f"resample-ref N={N} lam={lam} p={p:.4f}: mean {w.mean():.5f} var emp {w.var():.4e} formula {exr:.4e}")
# mean transient
N,lam,p=2000,0.05,0.03
w=np.full(40000,0.03*0+0.0); a=1-lam; m3=3*0+int(round(p*lam*N))
w[:]=0.5*p # E[w0]=0.015
m0=w[0]
for t in range(1,31):
    w=(rng.binomial(int(a*N),w)+m3)/N
print(f"transient t=30: emp mean {w.mean():.6f} formula {p+a**30*(m0-p):.6f}")

print("\n=== R4-B  lineage-resolved population operator on grid ===")
dx=0.01; x=np.arange(-30,45,dx)
mus=[-3,2,8]; sds=[1,0.8,0.6]; pis=[0.70,0.27,0.03]
def dens(mu,s): return norm.pdf(x,mu,s)
comp=[pis[k]*dens(mus[k],sds[k]) for k in range(3)]
b12,b23=-0.10488,5.58249
reg=[x<b12,(x>=b12)&(x<b23),x>=b23]
def mass(f,mask): return f[mask].sum()*dx
def run(lam,h,sig,T,refcomp):
    a=1-lam; s=np.sqrt(h*h+sig*sig)
    # q0 = P0 comp * N_h
    q=[gaussian_filter1d(c,h/dx,mode='constant') for c in comp]
    Rh=[gaussian_filter1d(c,h/dx,mode='constant') for c in refcomp]
    for t in range(T):
        q=[gaussian_filter1d(a*q[k],s/dx,mode='constant')+lam*Rh[k] for k in range(3)]
    return q
def analytic_op(lam,h,sig,T,refcomp):
    # same recursion on total, independent code path
    a=1-lam; s=np.sqrt(h*h+sig*sig)
    Q=gaussian_filter1d(sum(comp),h/dx,mode='constant'); R=gaussian_filter1d(sum(refcomp),h/dx,mode='constant')
    for t in range(T): Q=gaussian_filter1d(a*Q,s/dx,mode='constant')+lam*R
    return Q
for lam,h,T in [(0.05,0.6,59),(0.02,0.3,149)]:
    for label,rc in [("target-averaged",comp)]:
        q=run(lam,h,0.0,T,rc); Q=analytic_op(lam,h,0.0,T,rc)
        S=sum(q); print(f"lam={lam} h={h} T={T} {label}: max|sum_k q_k - Q| = {np.abs(S-Q).max():.2e}")
        A=np.array([[mass(q[k],reg[r]) for r in range(3)] for k in range(3)])
        print(" A[k,r] (rows=ancestry, cols=region):"); print(np.round(A,5))
        print(" row sums",np.round(A.sum(1),5)," col sums",np.round(A.sum(0),5)," region mass of Q",np.round([mass(Q,reg[r]) for r in range(3)],5))
        P0reg=[mass(sum(comp),reg[r]) for r in range(3)]
        print(" region mass of P0:",np.round(P0reg,5))
# seed-conditional: a reference draw of m points, labelled
def ref_draw(m,rs):
    k=rs.choice(3,size=m,p=pis); r=rs.normal(np.array(mus)[k],np.array(sds)[k]); return k,r
def refcomp_from(k,r,m):
    out=[]
    for j in range(3):
        pts=r[k==j]; d=np.zeros_like(x)
        idx=np.round((pts-x[0])/dx).astype(int); np.add.at(d,idx,1.0/(m*dx))  # atoms as grid spikes (density units)
        out.append(d)
    return out
lam,h,m=0.05,0.6,100; T=59; a=1-lam
rs=np.random.default_rng(7)
for seed in range(3):
    k,r=ref_draw(m,rs); ph=[(k==j).mean() for j in range(3)]
    rc=refcomp_from(k,r,m)
    q=run(lam,h,0.0,T,rc)
    mk=[q[j].sum()*dx for j in range(3)]
    pred=[ph[j]+a**T*(pis[j]-ph[j]) for j in range(3)]
    print(f"seed {seed}: pi_hat={np.round(ph,3)}  mass q^(k)_T={np.round(mk,5)}  predicted pi_hat+a^T(pi-pi_hat)={np.round(pred,5)}")
    Q=analytic_op(lam,h,0.0,T,rc); print("   max|sum q - Q_cond| =",f"{np.abs(sum(q)-Q).max():.2e}")

print("\n=== R4-C  seed-conditional fixed point vs target-averaged: W1 bound ===")
def cdf_from_density(d): return np.cumsum(d)*dx
def W1(d1,d2): return np.abs(cdf_from_density(d1)-cdf_from_density(d2)).sum()*dx
F0=sum(comp)
for lam,h in [(0.05,0.6),(0.02,0.3)]:
    a=1-lam; T=int(np.ceil(np.log(1e-7)/np.log(a))) if lam>0.04 else 800
    qs=run(lam,h,0.0,T,comp); Qs=sum(qs)
    rs=np.random.default_rng(11); rows=[]
    for m in [20,50,100,200,400]:
        d=[];dq=[]
        for rep in range(40):
            k,r=ref_draw(m,rs); rc=refcomp_from(k,r,m)
            Fr=sum(rc)
            w_ref=W1(Fr,F0)
            Qr=analytic_op(lam,h,0.0,T if lam>0.04 else 800,rc)
            d.append(w_ref); dq.append(W1(Qr,Qs))
        d=np.array(d);dq=np.array(dq)
        rows.append((m,d.mean(),dq.mean(),np.all(dq<=d+1e-3)))
    print(f"lam={lam} h={h}: m, E W1(Pref,P0), E W1(Q*_ref,Q*), bound holds in all draws")
    for r_ in rows: print("  ",r_[0],f"{r_[1]:.4f} {r_[2]:.4f}",r_[3])
