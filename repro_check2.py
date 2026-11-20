import pickle,time
import numpy as np
import phase3_pilot as pp
from repro_check import cmp,strip
saved={r['seed']:r for r in pickle.load(open('phase3_pilot_stage1.pkl','rb'))}[1]
t=time.time(); new=pp.job(1); mx,bad=cmp(strip(new),strip(saved)); print(f'stage1 seed 1: max abs diff {mx:.3e} issues {bad[:3]} ({time.time()-t:.0f}s)',flush=True)
