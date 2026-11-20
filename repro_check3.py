import pickle,time,numpy as np
import gmm_phase2d as d
from repro_check import cmp
R=pickle.load(open('gmm_phase2d_results.pkl','rb')); saved={(r[0],r[1]):r[2] for r in R}[(100,0)]
t=time.time(); new=d.job((100,0))[2]; mx,bad=cmp(new,saved); print(f'gmm2d (m=100,seed 0): max abs diff {mx:.3e} issues {bad[:3]} ({time.time()-t:.0f}s)',flush=True)
