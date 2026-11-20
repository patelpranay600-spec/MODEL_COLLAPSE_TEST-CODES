import pickle, numpy as np
from scipy import stats
R=sorted(pickle.load(open('phase3_s27.pkl','rb')),key=lambda r:r['seed']); n=len(R); print('n',n)
def v(key,win,q): return np.array([r['cells'][key][win][q] for r in R])
def ms(a): return f'{a.mean():8.4f}±{a.std(ddof=1)/np.sqrt(len(a)):.4f}'
def holm(ps):
    o=np.argsort(ps); adj=np.empty(len(ps)); run=0
    for rank,i in enumerate(o): run=max(run,(len(ps)-rank)*ps[i]); adj[i]=min(1,run)
    return adj
E=[('W1','W1'),('etau','e_tau'),('eR','e_R')]
fam={'F1 (m at lam=.05)':[(m,0.05) for m in (20,50,100,200,400)],'F2 (lam at m=100)':[(100,l) for l in (0.01,0.02,0.05,0.10)],'F3 (diagonal)':[(20,0.01),(40,0.02),(100,0.05),(200,0.10)]}
for win in ['R','G']:
    print(f'\n===== window {win} ({"each cell T5 window" if win=="R" else "generations 16-20"}) =====')
    for fn,cells in fam.items():
        print(f'\n{fn}')
        print(' (m,lam) T5 | W1 | W1_refonly(R_h) | W1(Qref,P0) | W1(C,Qref) | e_tau | e_R | |eps_tot| | pihat==0 seeds')
        for c in cells:
            et=np.abs((np.array([r['cells'][c]['om3_T5'] for r in R])-0.03)/0.03) if win=='R' else None
            ph=np.array([r['cells'][c]['pihat'] for r in R]); eg=np.abs((v(c,win,'om3')-0.03)/0.03)
            print(f' {str(c):12s} {R[0]["cells"][c]["T"]:3d} | {ms(v(c,win,"W1"))} | {ms(np.array([r["cells"][c]["Rh"]["W1"] for r in R]))} | {ms(np.array([r["cells"][c]["W1_Qref_P0"] for r in R]))} | {ms(np.array([r["cells"][c]["W1_C_Qref"] for r in R]))} | {ms(v(c,win,"etau"))} | {ms(v(c,win,"eR"))} | {ms(eg)} | {(ph==0).sum()}')
        print(' adjacent paired contrasts (later - earlier), Holm within endpoint x family:')
        for q,nm in E+[('eps','|eps_tot|')]:
            ps=[];ds=[]
            for a,b in zip(cells[:-1],cells[1:]):
                if q=='eps': d=np.abs((v(b,win,'om3')-0.03)/0.03)-np.abs((v(a,win,'om3')-0.03)/0.03)
                else: d=v(b,win,q)-v(a,win,q)
                ds.append(d); ps.append(stats.ttest_1samp(d,0)[1])
            adj=holm(np.array(ps)); print(f'  {nm:9s}: '+'; '.join(f'{a}->{b}: {d.mean():+.4f} (Holm p={p:.3f})' for (a,b),d,p in zip(zip(cells[:-1],cells[1:]),ds,adj)))
print('\nProp 3 check at T5 (conditional on retained set): mean (omega_T5 - pihat)^2 vs mean V^wr(pihat) = pihat(1-pihat)/(N lam (2-lam)+(1-lam))')
for c in [(20,.05),(50,.05),(100,.05),(200,.05),(400,.05),(100,.01),(100,.02),(100,.10)]:
    ph=np.array([r['cells'][c]['pihat'] for r in R]); om=np.array([r['cells'][c]['om3_T5'] for r in R]); lam=c[1]; vw=ph*(1-ph)/(2000*lam*(2-lam)+(1-lam)); print(f' {c}: observed {np.mean((om-ph)**2):.3e}  predicted {vw.mean():.3e}  ratio {np.mean((om-ph)**2)/vw.mean():.2f}')
