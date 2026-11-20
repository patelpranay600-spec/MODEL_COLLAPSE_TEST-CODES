import pickle,numpy as np,sys,gmm_phase2d as d,gmm_phase2 as g
from multiprocessing import Pool
VI=d.VARS.index('V3')
def j(a):
    m,seed=a; rng=np.random.default_rng(d.BASE2+seed*1000+m); k=rng.choice(3,size=m,p=g.PI); ref=rng.normal(g.MU[k],g.SD[k]); rs=(d.BASE2+seed)*7%1000000000; mod0=d.fitv(ref,'V0',rs,None); out={}
    for reg in ['R2','R3']:
        r=np.random.default_rng(d.BASE2+500000+seed*100+m+VI*7+(reg=='R3')*3); mod=mod0; tr=[]
        for t in range(1,d.T+1):
            D=g.sample(mod,d.NS if reg=='R3' else d.N,r)
            if reg=='R3': D=np.concatenate([D,ref[r.integers(0,m,d.NR)]])
            mod=d.fitv(D,'V3',(rs+t)%4294967295,mod); tr.append(g.metrics(mod))
        out[reg]=np.array(tr)
    return seed,out
if __name__=="__main__":
    R=pickle.load(open('gmm_phase2d_results.pkl','rb')); S={(r[0],r[1]):r[2] for r in R}; iW=g.NAMES.index('W1'); seeds=list(range(int(sys.argv[1])))
    with Pool(2) as p: res=p.map(j,[(100,s) for s in seeds])
    for reg in ['R2','R3']:
        a=np.array([S[(100,s)][('V3',reg)][-1,iW] for s in seeds]); b=np.array([o[reg][-1,iW] for s,o in res]); same=(np.array([S[(100,s)][('V3',reg)] for s in seeds])==np.array([o[reg] for s,o in res])).all((1,2))
        print(f'V3 {reg} m=100 seeds0-{len(seeds)-1}: identical trajectories {same.sum()}/{len(seeds)}; mean W1@T saved {a.mean():.3f} new {b.mean():.3f}; mean paired |diff| {np.abs(a-b).mean():.3f}; corr {np.corrcoef(a,b)[0,1]:.2f}',flush=True)
