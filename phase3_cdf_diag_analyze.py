import pickle,numpy as np
D=pickle.load(open('phase3_cdf_diag.pkl','rb')); D=sorted(D,key=lambda r:r['seed']); n=len(D); S1={r['seed']:r for r in pickle.load(open('phase3_pilot_stage1.pkl','rb'))}
def ci(x): x=np.asarray(x,float); return x.mean(),1.96*x.std(ddof=1)/np.sqrt(len(x))
rw=lambda r,arm,key: np.mean([r[arm][t-1][key] for t in range(55,60)],0)
print(f"{n} seeds. DIAGNOSTIC on previously used pilot seeds; not confirmatory.")
print("\n=== reproducibility vs stage 1 (R-window mean W1, tau, rv; max abs difference over seeds) ===")
dW=[];dT=[];dR=[];cW=[]
for r in D:
    s=S1[r['seed']]; key=('1.0',0.2,20,20) if ('1.0',0.2,20,20) in s['B'] else ('smin',0.2,20,20)
    dW.append(abs(rw(r,'B','W1')-s['B'][key]['W1'])); dT.append(abs(rw(r,'B','tau')-s['B'][key]['tau'])); dR.append(abs(rw(r,'B','rv')-s['B'][key]['rv'])); cW.append(abs(rw(r,'C','W1')-s['C']['W1']))
print(f" B: max|dW1| {max(dW):.2e}, max|dtau| {max(dT):.2e}, max|drv| {max(dR):.2e} ; C: max|dW1| {max(cW):.2e}")
print("\n=== regional W1 contributions (mean over seeds ± 95% HW): R1 x<-0.105 | R2 | R3 x>=5.58 | total | sum-of-regions check ===")
def row(lbl,v): 
    m=np.array([ci([x[k] for x in v]) for k in range(3)]); print(f" {lbl:22s} "+" | ".join(f"{m[k,0]:.4f}±{m[k,1]:.4f}" for k in range(3))+f" | total {np.mean([x.sum() for x in v]):.4f}")
for t in [1,5,10,20,40,59]:
    print(f"-- generation {t}")
    row("C",[r['C'][t-1]['reg'] for r in D]); row("B (1.0, .2, 20x20)",[r['B'][t-1]['reg'] for r in D])
print("-- R window (55-59 mean)"); row("C",[rw(r,'C','reg') for r in D]); row("B",[rw(r,'B','reg') for r in D])
for nm in ['R_h','R_match','R_loo']: row(nm+" (reference-only)",[r[nm]['reg'] for r in D])
print("total vs sum of regions, max abs error over all stored (law, generation):",max(abs(x['W1']-x['reg'].sum()) for r in D for x in r['B']+r['C']))
print("\n=== paired differences in regional W1 (R window; mean ± 95% HW) ===")
for a,b,lab in [('B','R_h','B - R_h'),('B','R_match','B - R_match'),('B','C','B - C')]:
    fa=(lambda r:rw(r,'B','reg')) ; 
    fb=(lambda r:r[b]['reg']) if b!='C' else (lambda r:rw(r,'C','reg'))
    d=np.array([fa(r)-fb(r) for r in D]); m=[ci(d[:,k]) for k in range(3)]; t=ci(d.sum(1)); print(f" {lab:13s} R1 {m[0][0]:+.4f}±{m[0][1]:.4f} | R2 {m[1][0]:+.4f}±{m[1][1]:.4f} | R3 {m[2][0]:+.4f}±{m[2][1]:.4f} | total {t[0]:+.4f}±{t[1]:.4f}")
print("\n=== generation-wise paired B - R_h and B - C, total W1 ===")
for t in [1,5,10,20,40,59]:
    a=ci([r['B'][t-1]['W1']-r['R_h']['W1'] for r in D]); b=ci([r['B'][t-1]['W1']-r['C'][t-1]['W1'] for r in D]); print(f" t={t:2d}: B-R_h {a[0]:+.3f}±{a[1]:.3f} ; B-C {b[0]:+.3f}±{b[1]:.3f} ; C W1 {np.mean([r['C'][t-1]['W1'] for r in D]):.3f} B W1 {np.mean([r['B'][t-1]['W1'] for r in D]):.3f} R_h {np.mean([r['R_h']['W1'] for r in D]):.3f}")
print("\n=== signed contribution ∫(F-F0)dx by region, R window (positive = law's CDF above the target: mass shifted left of the target there) ===")
for lbl,f in [('C',lambda r:rw(r,'C','sreg')),('B',lambda r:rw(r,'B','sreg')),('R_h',lambda r:r['R_h']['sreg']),('R_match',lambda r:r['R_match']['sreg'])]:
    v=np.array([f(r) for r in D]); print(f" {lbl:8s} "+" | ".join(f"{v[:,k].mean():+.4f}" for k in range(3)))
print("\n=== matched-smoothing hypothesis (post hoc): B vs reference-only KDE at sqrt(h^2+s_min^2) (R_match), R window ===")
for key,lab in [('tau','tail mass'),('rv','variance ratio'),('rho3','region-3 mass'),('W1','W1')]:
    B=np.array([rw(r,'B',key) for r in D]); M=np.array([r['R_match'][key] for r in D]); H_=np.array([r['R_h'][key] for r in D]); C=np.array([rw(r,'C',key) for r in D])
    d=ci(B-M); print(f" {lab:15s} B {B.mean():.4f} | R_match {M.mean():.4f} | R_h {H_.mean():.4f} | C {C.mean():.4f} | B-R_match {d[0]:+.4f}±{d[1]:.4f} | seed corr(B,R_match) {np.corrcoef(B,M)[0,1]:+.2f}")
print(" share of B's excess tail over R_h explained by R_match: (R_match-R_h)/(B-R_h) = %.2f"%((np.mean([r['R_match']['tau'] for r in D])-np.mean([r['R_h']['tau'] for r in D]))/(np.mean([rw(r,'B','tau') for r in D])-np.mean([r['R_h']['tau'] for r in D]))))
print(" share of B's W1 gap to R_h explained by R_match: (R_match-R_h)/(B-R_h) = %.2f"%((np.mean([r['R_match']['W1'] for r in D])-np.mean([r['R_h']['W1'] for r in D]))/(np.mean([rw(r,'B','W1') for r in D])-np.mean([r['R_h']['W1'] for r in D]))))
