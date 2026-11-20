import pickle,sys,numpy as np
tag=sys.argv[1]; R=pickle.load(open(f'phase3_pilot_{tag}.pkl','rb')); n=len(R)
def wilson(k,n,z=1.96):
    if n==0: return (np.nan,np.nan)
    p=k/n; d=1+z*z/n; c=(p+z*z/(2*n))/d; h=z*np.sqrt(p*(1-p)/n+z*z/(4*n*n))/d; return c-h,c+h
cfgs=sorted({k for r in R for k in r['B']},key=lambda k:(k[2],k[1],['smin','1.0','2.5'].index(k[0])))
print(f"{tag}: {n} seeds. s_min median {np.median([r['smin'] for r in R]):.3f} range {min(r['smin'] for r in R):.2f}-{max(r['smin'] for r in R):.2f}; refs with no rare point: {sum(not r['has_rare'] for r in R)}")
print("rule changed a nominal s0 (s_min > nominal):",{v:sum(any(x[0]==v for x in r['nominal_changed']) for r in R) for v in ['1.0','2.5']}, "; deduplicated settings (seeds with <3 s0 values):",sum(len(r['s0'])<3 for r in R))
print("tail mass reference points: target 0.004533 | smoothing floor P0*N_0.6 0.01238 | amended ceiling 0.0248 | R_h mean %.4f | R_loo mean %.4f"%(np.mean([r["R_h"]["tau"] for r in R]),np.mean([r["R_loo"]["tau"] for r in R])))
print("C (uncorrected) R-window means: W1 %.3f rho3 %.4f tau %.4f rv %.3f | R_h W1 %.3f | R_loo W1 %.3f"%tuple(np.mean(x) for x in ([r['C']['W1'] for r in R],[r['C']['rho3'] for r in R],[r['C']['tau'] for r in R],[r['C']['rv'] for r in R],[r['R_h']['W1'] for r in R],[r['R_loo']['W1'] for r in R])))
def J(o):
    A=o['A'];w=o['w']; 
    if np.isnan(A[2]).any(): return np.nan,np.nan,np.nan
    inflow=sum(w[i]*A[i,2] for i in (0,1) if not np.isnan(A[i,2])); return inflow-w[2]*(1-A[2,2]),inflow,A[2,2]
rows=[]
print("\nGUARDRAIL failure rates (fraction of evaluated seeds; Wilson 95%) — variance | tail | region-3 | any ; tail ceiling AMENDED to 0.0248 (Section 20.B; original 0.009 shown as 'orig tail'); elimination if any single lower bound > 0.10")
for c in cfgs:
    S=[r for r in R if c in r['B']]; fv=ft=f3=fa=ft0=0; fin=True
    for r in S:
        b=r['B'][c]; C=r['C']; v=abs(b['rv']-1)>abs(C['rv']-1)+0.05; t=b['tau']>0.0248; t0=b['tau']>0.009; ft0+=t0; q=not(0.015<=b['rho3']<=0.060); fv+=v; ft+=t; f3+=q; fa+=(v or t or q); fin&=b['finite']
    k=len(S); lo=[wilson(x,k)[0] for x in (fv,ft,f3)]; elim=any(l>0.10 for l in lo)
    rows.append((c,k,fv,ft,f3,fa,elim,fin))
    print(f" s0={c[0]:4s} c={c[1]:<4} L,K={c[2]},{c[3]:<2} n={k:2d} | var {fv/k:.2f} [{wilson(fv,k)[0]:.2f},{wilson(fv,k)[1]:.2f}] | tail {ft/k:.2f} [{wilson(ft,k)[0]:.2f},{wilson(ft,k)[1]:.2f}] | R3 {f3/k:.2f} [{wilson(f3,k)[0]:.2f},{wilson(f3,k)[1]:.2f}] | any {fa/k:.2f} | orig tail(0.009) {ft0/k:.2f} | finite {fin} | {'ELIMINATE' if elim else 'keep'}")
print("\nENDPOINTS (R-window seed means) and paired differences vs C and R_h (B-C, B-R_h; mean ± se), rare-lineage transport")
print(" config                     | W1    | B-C            | B-R_h          | rho3   | tau    | rv    | A33 B / C  | inflow B / C | J3 B / C | wall s | updates")
for c in cfgs:
    S=[r for r in R if c in r['B']]; b=lambda k:[r['B'][c][k] for r in S]; ms=lambda x:(np.mean(x),np.std(x,ddof=1)/np.sqrt(len(x)) if len(x)>1 else np.nan)
    dC=ms(np.array(b('W1'))-np.array([r['C']['W1'] for r in S])); dR=ms(np.array(b('W1'))-np.array([r['R_h']['W1'] for r in S]))
    tb=[J(r['B'][c]) for r in S]; tc=[J(r['C']) for r in S]; g=lambda L,i: np.nanmean([x[i] for x in L]) if not np.all(np.isnan([x[i] for x in L])) else np.nan
    print(f" s0={c[0]:4s} c={c[1]:<4} LK={c[2]:>2}x{c[3]:<2} | {np.mean(b('W1')):.3f} | {dC[0]:+.3f}±{dC[1]:.3f} | {dR[0]:+.3f}±{dR[1]:.3f} | {np.mean(b('rho3')):.4f} | {np.mean(b('tau')):.4f} | {np.mean(b('rv')):.3f} | {g(tb,2):.2f}/{g(tc,2):.2f} | {g(tb,1):.3f}/{g(tc,1):.3f} | {g(tb,0):+.3f}/{g(tc,0):+.3f} | {np.mean(b('time')):5.0f} | {c[2]*c[3]}")
print("(A33, inflow, J3 averaged over seeds whose reference/lineage composition makes them defined; reference-set SD not shown here)")
