import json
db={
"meta":{"project":"Model collapse in recursive synthetic data: noise + reference-guided correction","checkpoint_date":"2026-11-10","protocol_file":"model_collapse_math_protocol.md (Revision 4 + addenda, through Sections 25.B and DRAFT 26-28)","evidence_tags":{"[A]":"analytic","[G]":"grid/numeric exact","[MC]":"Monte Carlo","[H]":"heuristic","[U]":"untested/unvalidated"}},
"setting":{"target":"0.70 N(-3,1) + 0.27 N(2,0.8^2) + 0.03 N(8,0.6^2)","recursion":"D_t=(1-lam)*synthetic+lam*reference, replace-style","lam":0.05,"N":2000,"n_ref":100,"reference_draw":"with replacement from retained set m=100","h":0.6,"T5":59,"R_window":"55-59","G_window":"16-20","regions":{"b12":-0.10488,"b23":5.58249},"target_tail_mass":0.004533,"target_region3_mass":0.0300,"smoothing_floor_tail_P0*N_0.6":0.01238},
"seed_registry":[
{"range":"1-30","use":"stage-1 pilot and CDF diagnostic (retired as diagnostic)"},
{"range":"31-60","use":"Block 1 (Section 26, run)"},
{"range":"61-90","use":"Block 2 factorial (Section 27, run)"},
{"range":"91-100","use":"reserved pilot buffer"},
{"range":"201-230","use":"Section 25 smoothing-control pilot (run, closed)"},
{"range":"231-260","use":"reserved for deferred Block 4 mixing diagnostic"},
{"range":"301-330","use":"Block 3 H1 feasibility (Section 28, run)"},
{"range":"1001-1300","use":"confirmatory H1 run (Section 30); NOT used for H2/AB"}],
"findings":[
{"id":"F1","phase":"1 (KDE)","claim":"Population recursion is a reference-only scale-mixture estimator; finite-population sampling adds major error; rare-mode failure driven by missing reference observations; recursive KDE does not beat a properly selected reference-only estimator.","tag":"[A]/[MC]","status":"frozen, successful negative result"},
{"id":"F2","phase":"2 (GMM)","claim":"Pure recursion degrades W1 in 9/9 cells; warm start helps; multi-start inconclusive; no variant removes the degradation (Phase 2d EM comparison, +-0.05 margin for 2d only).","tag":"[MC]","status":"frozen"},
{"id":"F3","phase":"3 validation","claim":"S-NP score exact (analytic vs 40-digit mpmath agreement 4e-16); Langevin stationary variance within 0.71%; lambda=0 variance recursion matches finite-N correction (1-1/N); all 8 arm invariants pass. V18(a) criterion not met as written (2e-5 vs 1e-6), kept labelled.","tag":"[A]/[G]/[MC]","status":"closed"},
{"id":"F4","phase":"3 pilot stage 1","claim":"0/12 S-NP Langevin schedule configs passed guardrails on seeds 1-30; region-3 absolute guard confounded by reference composition; B tail about 2x reference-only. Negative result scoped to the tested schedule family, not a rejection of score correction.","tag":"[MC]","status":"closed (Section 21.0)"},
{"id":"F5","phase":"3 mechanism","claim":"B approx reference-only KDE smoothed at sqrt(h^2+s_min^2) (R_match): B-R_match tail +0.0001+-0.0004 (corr 1.00), variance ratio +0.029, W1 +0.018 (seeds 1-30, s0>=1.0 regime).","tag":"[MC] post hoc","status":"intervention-tested in F6"},
{"id":"F6","phase":"Section 25 (seeds 201-230)","claim":"Matching total smoothing (learner h_fit=sqrt(max(0.15^2,h^2-s_min^2))) removes tail and variance excess on unclamped seeds: tail ratio 0.98 (100% pass; B_orig 1.52, 43% pass), rv 1.008. R_match explains 1.06 (CI 1.02-1.12) of B_orig tail excess over R_h.","tag":"[MC]","status":"confirmed for tail/variance"},
{"id":"F7","phase":"Section 25 (seeds 201-230)","claim":"B_match does NOT tie reference-only in W1: vs R_eff +0.027 (90% CI 0.006-0.049, TOST margin 0.040 fails), vs R_h +0.042 (0.016-0.068), vs R_loo +0.038. rho3 guard passes on 78% of unclamped seeds (<90%). 7 of 30 seeds clamped; tail guard fails there (14% pass).","tag":"[MC]","status":"refuted (tie bar not met)"},
{"id":"F8","phase":"Section 25 transient","claim":"W1 at t=1/5/10/20: C 0.260/0.523/0.684/0.851; B_orig 0.294/0.351/0.378/0.405; B_match 0.246/0.309/0.343/0.357; R_h 0.312. B_match plateaus near its R-window value 0.355 by t=20. The claim that 400 steps erase S_0 at t=1 is unsupported.","tag":"[MC]","status":"recorded"},
{"id":"F9","phase":"Section 25 residual","claim":"Cause of residual W1 gap unexplained. Hypotheses only: incomplete cross-region mixing; ULA discretisation bias (V21 variance bias 0.71%, too small on its face); h_fit<h learner effect.","tag":"[U]","status":"deferred (Block 4), needs fresh seeds 231-260 and frozen criteria"},
{"id":"F10","phase":"Section 25 vs C","claim":"B_orig - C = -0.647 (same learner); B_match - C = -0.678 (different learner, not a same-learner comparison). A33: C 0.85, B_orig 0.73, B_match 0.78 (reported only).","tag":"[MC]","status":"recorded"}],
"net_conclusion":"In the 1-D fixed-bandwidth KDE recursion with retained reference (m=100, lam=0.05), the tested reference-guided S-NP Langevin correction removes most uncorrected degradation but does not match or beat the reference-only benchmark; its tail/variance inflation is explained by double smoothing and removable by learner-side bandwidth matching on unclamped seeds. Scope: tested constant-level s_min Langevin; not score correction in general. Reference-only KDE remains the benchmark to beat (a benchmark, not a proven ceiling).",
"decisions_and_amendments":[
"Stabilisation: sustained 5%-of-gap with absolute floors (W1 0.001, rare-region mass 0.0001, variance ratio 0.001, log-tail 0.01), 10-generation window.",
"T_fix=60 for lambda=0; T_R=T5 primary; stabilisation windows sensitivity only.",
"+-0.05 W1 margin is for Phase 2d only; 0.040 TOST margin is pilot-specific (Section 25.A) and tentatively reused in Block 1 pending ratification.",
"Tail guard amendment (post-execution): 0.009 -> 0.0248 (2 x 0.0124) for h=0.6 only; audit trail kept (first 2 seeds under original).",
"C3 no-transport gate dropped; transition matrix is primary mechanistic read-out; s0=max(planned,s_min^LOO); independent unit = seed/reference set.",
"Section 25.A frozen and ratified (h_min=0.15, tail ratio<=1.35 on unclamped seeds, TOST 0.040, success bar 'tie with best reference-only benchmark'); closed in 25.B; no third S-NP parameter sweep.",
"Block 4 (mixing diagnostic) deferred."],
"audit_log":[
"Phase 2 write-up: two prose claims corrected against outputs.",
"Seed plan misstatement corrected (pilot 1-100, test 1001-1300 exist).",
"Pilot runs lost to sandbox restart twice; pilot script made resumable.",
"User's paired-W1 s0=1.0 'improvement' used different seed sets; paired data +0.007+-0.026 (no benefit).",
"Section 25 spec frozen and launched before ratification (ordering error; thresholds fixed before results seen); run interrupted at 14/30 and resumed; bootstrap RNG reseeding bug fixed (no threshold change); 3.29*0.067/sqrt(30)=0.0402 vs margin 0.040 noted.",
"Block 1 timing check: V acceptance diagnostic bug (always 0.5) fixed before any real run."],
"pending":[
"Paper compilation (user-side) from protocol Sections 21-30; test seeds 1001-1300 now used for H1 only.",
"Step 7 feasibility (paired SDs, power for P-H2a/b/c, Holm, variance guard for 300 test seeds).",
"Unvalidated: S-NN, DR calibration, DR* weight variance, equal-compute variants, peak memory.",
"Item 31 (P-H2c as fraction of C-to-oracle-R gap); optional lineage-transport test."],
"standing_constraints":[
"Do not launch test seeds 1001-1300 before schedule and analysis rules are locked; no retuning on test seeds.",
"Keep stability, distributional quality and ancestry separate.",
"Do not silently change frozen thresholds; record post-execution amendments with audit trail.",
"New questions need fresh seeds in a new namespace and pre-frozen criteria.",
"No new correction methods mid-experiment; no S-NP as a separate experiment.",
"Reference-only KDE is the benchmark to beat."],
"findings_blocks":[
{"id":"F11","phase":"Section 26.B (seeds 31-60)","claim":"No non-oracle baseline ties/beats matched reference-only in W1: DR-C -0.469, V-C -0.265, S-C -0.017 (ns); DR-R_h +0.196, V-R_h +0.400, S-R_h^S +0.802. R_h^S-R_h -0.154. DR* W1 0.156 < R_h. V and V* amputate rare mode (rho3 0.0056, 0.0039). DR ESS 9412/9500.","tag":"[MC]","status":"recorded; mechanisms [H]"},
{"id":"F12","phase":"Section 27.B (seeds 61-90)","claim":"W1 of C is set by fixed-point smoothing and by lambda (2.88, 1.74, 0.98, 0.66 for lam .01-.10), not by m for m>=50 (1.00 vs R_h 0.25->0.20). Prop 3 variance check approximate (ratios 0.45-1.56, n=30).","tag":"[MC]","status":"recorded"},
{"id":"F13","phase":"Section 28.B (seeds 301-330)","claim":"Coupled check passes (1.2e-14). At n=300, l_tau and W1 feasible for all kappa; eps_tot feasible at kappa 0.6,0.9 only; rho_v and eR secondary not feasible. Paired/unpaired SD ratio 0.61-1.04 (pairing modest).","tag":"[MC]","status":"feasibility only; H1 untested"}],
"findings_h1":[{"id":"F14","phase":"Section 29/30 (test seeds 1001-1300, H1 only)","claim":"H2/AB closed at pilot gate (no surviving candidate; no test-seed H2 contrast computed). Confirmatory H1: coupled check 1.4e-14; l_tau equivalent at all kappa; eps_tot equivalent at kappa 0.9 but not 0.6 (90% CI [0.006,0.102] vs margin 0.10); overall H1 NOT supported under the frozen IU rule, marginal and likely power (test s_d 0.50 vs pilot 0.42-0.45).","tag":"[MC]","status":"final; margins not widened"}],
"files":{
"protocol":"model_collapse_math_protocol.md; backups rev5_pre_*_backup.md",
"phase1_2":"research_record_phases1_2.md; gmm_phase2*.py/.md",
"phase3_validation":"phase3_val_a.py, phase3_val_a2.py, phase3_val_b*.py, phase3_step2_lam0.py, phase3_arms.py",
"pilot":"phase3_pilot.py, phase3_pilot_analyze.py, phase3_stage1_mech.py, phase3_cdf_diag*.py/.pkl",
"section25":"phase3_s25_pilot.py/.pkl, phase3_s25_analyze.py, phase3_s25_analysis.txt",
"blocks":"phase3_s26*.py/.pkl (done), phase3_s27*.py/.pkl (done), phase3_s28*.py/.pkl (done); outputs phase3_s2[678]_analysis.txt"}}
json.dump(db,open('model_collapse_checkpoint.json','w'),indent=1)
L=[f"# Model collapse research — context checkpoint ({db['meta']['checkpoint_date']})","",f"Protocol: {db['meta']['protocol_file']}. Evidence tags: "+"; ".join(f"{k} {v}" for k,v in db['meta']['evidence_tags'].items()),"","## Net conclusion",db['net_conclusion'],"","## Setting"]
L+=[f"- {k}: {v}" for k,v in db['setting'].items()]
L+=["","## Seed registry"]+[f"- {s['range']}: {s['use']}" for s in db['seed_registry']]
L+=["","## Findings"]+[f"- **{f['id']}** ({f['phase']}) {f['tag']} — {f['claim']} *Status: {f['status']}*" for f in db['findings']+db['findings_blocks']+db['findings_h1']]
for k,t in [('decisions_and_amendments','Decisions and amendments'),('audit_log','Audit log'),('pending','Pending'),('standing_constraints','Standing constraints')]: L+=["",f"## {t}"]+[f"- {x}" for x in db[k]]
L+=["","## File index"]+[f"- {k}: {v}" for k,v in db['files'].items()]
open('model_collapse_checkpoint.md','w').write("\n".join(L)+"\n")
