# Phase 2d — EM fitting-procedure comparison (pre-registered; written before any run)

Purpose: n_init > 1 (selection among local optima) and warm-starting (continuity of the optimisation path) target different failure mechanisms. Neither is assumed better. Phase 2c (exploratory, N = 20,000, 20 seeds) hinted warm start helps and tight EM alone does not; that is a hypothesis here, not a premise.

Seed namespace: BASE2 = 2026110900 (references, draws and EM seeds all new; no reuse of Phase 2/2b/2c seeds). Self-contained script `gmm_phase2d.py`.

Setting: target as before; learner K3 only (K2 never holds the rare component; BIC warm start is undefined when K changes). N = 2000, T = 60, 60 seeds, m ∈ {100, 400}. Regimes: R2 (λ = 0) primary; R3 (λ = 0.05 injection) secondary. Generation-0 fit on the reference uses the V0 procedure for all variants, so every variant starts from the same model.

Variants (generation t ≥ 1 fits):
- V0 default: cold k-means init, n_init = 1, tol 1e-3, max_iter 100 (Phase 2 procedure).
- V1 multi-start: cold, n_init = 5, tol 1e-6, max_iter 500; best log-likelihood kept.
- V2 warm: init from previous generation's parameters, n_init = 1, tol 1e-6, max_iter 500.
- V3 hybrid: best log-likelihood among the V2 warm fit and a V1 multi-start fit.

Endpoints (window = generations 51–60 means): W1, |ε_R|, rare-component presence (a fitted component with mean > b23 and weight ≥ 0.01), log tail ratio ℓ_τ (median). Primary endpoint: W1 growth = late − t=1 (paired by seed). Reference level: R0 = generation-0 fit W1.

Contrasts, paired by seed, 95% CI from paired SE (normal): V1−V0, V2−V0, V3−V0, V3−V1, V3−V2 on late W1 and on W1 growth; same under R3.

"≈" margin (judgement-based, flagged for ratification): |Δ late W1| ≤ 0.05 with the CI inside ±0.05.

Pre-specified interpretation:
| V1 vs V0 | V2 vs V0 | reading |
|---|---|---|
| better | ≈ | local-optimum selection is the main fitting failure |
| ≈ | better | path continuity is the main fitting failure |
| better | better | both contribute; V3 vs V1/V2 says whether they are additive |
| ≈ | ≈ | fitting procedure is not the driver; residual is finite-sample drift |
"better" = CI of late-W1 difference entirely below 0 (and below −0.05 to count as material).
Also reported regardless: whether any variant's R2 late W1 reaches within 0.05 of R0 (i.e. removes the degradation), and whether the variant effect survives under R3.

Not claimed: anything about N ≠ 2000, other T, other targets, BIC/K2 learners, or ancestry.
