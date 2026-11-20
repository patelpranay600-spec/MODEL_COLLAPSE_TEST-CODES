import pickle,numpy as np,gmm_phase2d as d,gmm_phase2 as g
R=pickle.load(open('gmm_phase2d_results.pkl','rb')); saved={(r[0],r[1]):r[2] for r in R}[(100,0)]; new=d.job((100,0))[2]
print('metrics',g.NAMES)
for k in saved:
    if k==('R0',): print(k,'diff',np.abs(np.array(list(saved[k].values())) -np.array(list(new[k].values()))).max() if isinstance(saved[k],dict) else np.abs(np.asarray(saved[k],float)-np.asarray(new[k],float)).max()); continue
    a,b=saved[k],new[k]; df=np.abs(a-b); first=int(np.argmax((df>1e-9).any(1))) if (df>1e-9).any() else None
    iW=g.NAMES.index('W1'); print(k,'first divergent gen',first,'max diff',round(float(df.max()),3),'W1 at T saved/new',round(float(a[-1,iW]),3),round(float(b[-1,iW]),3))
