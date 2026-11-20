# Re-runs one seed of each Phase-3 block from scratch and compares with the saved pickle (exact and 1e-12).
import pickle, numpy as np, sys, time, importlib
def cmp(a,b,path=''):
    mx=0.0; bad=[]
    if isinstance(a,dict):
        if set(map(str,a))!=set(map(str,b)): return 0.0,[path+' keys differ']
        for k in a:
            kb=k if k in b else next(x for x in b if str(x)==str(k)); m,bd=cmp(a[k],b[kb],path+'/'+str(k)); mx=max(mx,m); bad+=bd
    elif isinstance(a,(list,tuple)):
        if len(a)!=len(b): return 0.0,[path+' len']
        for i,(x,y) in enumerate(zip(a,b)): m,bd=cmp(x,y,path+f'[{i}]'); mx=max(mx,m); bad+=bd
    elif isinstance(a,np.ndarray) or isinstance(a,(float,np.floating,int,np.integer)):
        x=np.asarray(a,float); y=np.asarray(b,float)
        if x.shape!=y.shape: return 0.0,[path+' shape']
        d=np.nanmax(np.abs(x-y)) if x.size else 0.0
        if np.isnan(x).sum()!=np.isnan(y).sum(): bad.append(path+' nan pattern')
        mx=float(d)
    elif a is None or isinstance(a,(bool,str)):
        if a!=b: bad.append(path+' value')
    return mx,bad
def strip(o):
    if isinstance(o,dict): return {k:strip(v) for k,v in o.items() if k not in('time','cnt')}
    if isinstance(o,list): return [strip(x) for x in o]
    return o
jobs={'s25':('phase3_s25_pilot','phase3_s25_pilot.pkl',201),'s26':('phase3_s26','phase3_s26.pkl',31),'s27':('phase3_s27','phase3_s27.pkl',61),'s28':('phase3_s28','phase3_s28.pkl',301),'s30':('phase3_s28','phase3_s30.pkl',1001)}
for name in (sys.argv[1:] if __name__=="__main__" else []):
    mod,pk,seed=jobs[name]; m=importlib.import_module(mod); saved={r['seed']:r for r in pickle.load(open(pk,'rb'))}[seed]; t=time.time(); new=m.job(seed)
    mx,bad=cmp(strip(new),strip(saved)); print(f'{name} seed {seed}: max abs diff {mx:.3e}  structural issues {bad[:3]}  ({time.time()-t:.1f}s)',flush=True)
