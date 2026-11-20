# Pickle-independent per-seed exports of the Phase-3 blocks (plain CSV).
import pickle,csv,numpy as np
def w(fn,rows,cols):
    with open(fn,'w',newline='') as f: c=csv.writer(f); c.writerow(cols); c.writerows(rows)
A=sorted(pickle.load(open('phase3_s25_pilot.pkl','rb')),key=lambda r:r['seed']); rows=[]
for r in A:
    for a in ['C','R_h','R_loo','R_match','R_eff','B_orig','B_match']: x=r[a]; rows.append([r['seed'],a,r['smin'],int(r['clamped']),x['W1'],x['rho3'],x['tau'],x['rv']])
w('export_s25_per_seed.csv',rows,['seed','arm','smin','clamped','W1','rho3','tau','rv'])
B=sorted(pickle.load(open('phase3_s26.pkl','rb')),key=lambda r:r['seed']); rows=[]
for r in B:
    for a in ['C','S','DR','V','DR*','V*','R_h','R_hS']: x=r[a]; rows.append([r['seed'],a,x['W1'],x['rho3'],x['tau'],x['rv']])
w('export_s26_per_seed.csv',rows,['seed','arm','W1','rho3','tau','rv'])
C=sorted(pickle.load(open('phase3_s27.pkl','rb')),key=lambda r:r['seed']); rows=[]
for r in C:
    for (m,l),c in sorted(r['cells'].items()):
        for win in ['R','G']: x=c[win]; rows.append([r['seed'],m,l,win,x['W1'],x['rho3'],x['eR'],x['tau'],x['etau'],x['rv'],x['om3'],c['pihat'],c['Rh']['W1'],c['W1_Qref_P0'],c['W1_C_Qref']])
w('export_s27_per_seed.csv',rows,['seed','m','lam','window','W1','rho3','eR','tau','etau','rv','om3','pihat','W1_Rh','W1_Qref_P0','W1_C_Qref'])
for tag,fn in [('s28','phase3_s28.pkl'),('s30','phase3_s30.pkl')]:
    D=sorted(pickle.load(open(fn,'rb')),key=lambda r:r['seed']); rows=[]
    for r in D:
        for kap,v in r['k'].items():
            for arm in ['N','H']: x=v[arm]; rows.append([r['seed'],kap,arm,x['W1'],x['ltau'],x['rv'],x['eR'],x['etot'],v['coupled_maxdiff']])
    w(f'export_{tag}_per_seed.csv',rows,['seed','kappa','arm','W1','ltau','rv','eR','etot','coupled_maxdiff'])
print('ok')
