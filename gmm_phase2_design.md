# Phase 2 — controlled-capacity GMM recursion: pre-registered design (written before any run)

Seed namespace: 2026110800 + seed (new; disjoint from Phase 1). Target P0, N = 2000 per generation, T = 60 generations, 60 paired seeds, reference sizes m ∈ {20, 100, 400}. One retained reference set per (seed, m), shared across all learners and regimes.

**Learners (sklearn Gaussian mixture, EM, kmeans init, n_init = 1, max_iter = 100, reg_covar = 1e-6):**
- K3: correct component count (reported as the favourable case, not as a realistic learner).
- K2: under-specified capacity.
- BIC: K ∈ {1,…,6} chosen by BIC at every generation (data-driven, no target knowledge).

**Regimes:**
- R0 reference-only: fit on the m reference points.
- R2 pure recursion (λ = 0): D_t = N draws from model_{t−1}; model_0 = R0. R1 (one-generation synthetic) is R2 at t = 1.
- R3 injection (λ = 0.05): D_t = 1900 draws from model_{t−1} + 100 draws with replacement from the reference.
- R4 fresh-data control: fit on N iid draws from P0 (equal sample budget, ideal information).
- R5 reference-bootstrap control: fit on N draws with replacement from the reference (equal sample budget, reference information only, no generative model).

**Endpoints (per generation):** W1 to P0 (analytic mixture CDF on a grid); signed ε_R = (m_R − 0.03)/0.03 for the region above b23 = 5.58249 (analytic); ℓ_τ = log(τ/τ0) with τ the mass in (−∞,−5.5] ∪ [9.5,∞) (analytic; GMM tails are never exactly zero); variance ratio ρ_v; mode count (basins of the fitted density with mass ≥ 0.005); rare-component presence (a fitted component with mean above b23 and weight ≥ 0.01); selected K (BIC). **Ancestry is not tracked:** generated points carry no component label through a GMM fit.

**Primary window:** generations 51–60 (mean of per-generation metrics). Trajectories at t = 1, 5, 10, 20, 30, 40, 50, 60.

**Pre-specified contrasts (paired by seed, mean ± se):** (a) R2 late − R2 at t = 1 (growth over generations); (b) OLS slope of the metric on generation within R2 and R3; (c) R2 late − R5 and R2 late − R0 (is recursion worse than finite reference alone?); (d) R3 late − R2 late (effect of injection); (e) R2 late − R4.

**Exit condition for "genuine, reproducible degradation attributable to the learner":** for some learner and m, (a) and (b) on W1 or on rare-region mass show degradation with the paired 95% CI excluding 0 **and** (c) R2 late is worse than R5/R0 with the paired CI excluding 0. Results are reported for all learners and m, including those that do not meet it. No learner or setting is dropped because it fails to show degradation, and none is chosen for showing it.

**Limits fixed in advance:** one N, one T, λ = 0.05 only for injection; no correction method is involved; mode-count and rare-presence definitions are as above.
