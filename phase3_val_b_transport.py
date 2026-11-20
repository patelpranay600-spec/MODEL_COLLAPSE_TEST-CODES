# Full signed conditional transition matrices from the V22 results (no new run). Rows = input label (conditional, sum 1), cols = output region.
import pickle, numpy as np
R=pickle.load(open('phase3_val_b_results.pkl','rb'))
def cond(A): return A/A.sum(1,keepdims=True)
print("Conditional dA = A_B - A_in (rows: input label 1..3; cols: output region R1..R3). C comparator = no correction (output = input). Mean over 10 references; reference-level SD of J3 and A33 in brackets.")
for m in [20,100,400]:
    RR=[r for r in R if r['m']==m]
    for kind in ['Qstar']:
        for s0 in ['smin','1.0']:
            for c in [0.05,0.2]:
                dAs=[];J3=[];A33=[];A33C=[];inflow=[];outflow=[]
                for r in RR:
                    Ab=r[(kind,s0,c)]['A']; Ai=r[(kind,'in')]['A']; w=Ai.sum(1); cb=cond(Ab); ci=cond(Ai)
                    dAs.append(cb-ci); inflow.append(sum(w[i]*cb[i,2] for i in (0,1))); outflow.append(w[2]*(1-cb[2,2])); J3.append(inflow[-1]-outflow[-1]); A33.append(cb[2,2]); A33C.append(ci[2,2])
                d=np.mean(dAs,0)
                print(f"\nm={m} {kind} s0={s0} c={c}\n  dA=\n"+"\n".join("   "+" ".join(f"{x:+.3f}" for x in row) for row in d))
                print(f"  J3 {np.mean(J3):+.4f} [{np.std(J3,ddof=1):.4f}]  inflow(non-rare->R3) {np.mean(inflow):.4f}  outflow(rare leaving R3) {np.mean(outflow):.4f}  A33 {np.mean(A33):.3f} [{np.std(A33,ddof=1):.3f}] vs C {np.mean(A33C):.3f}")
