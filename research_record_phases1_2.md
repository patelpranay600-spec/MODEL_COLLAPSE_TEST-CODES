# Consolidated research record — Phases 1, 2, 2b/2c, 2d and protocol closure

All conclusions are scoped to the tested target (1-D three-component Gaussian mixture, weights 0.70/0.27/0.03), learners (Gaussian KDE; sklearn EM-GMM K3/K2/BIC), N, and horizons stated in each part. No correction method (S-NP or reference-guided score) has been run.

| Phase | Status |
|---|---|
| 1 KDE baseline | Frozen |
| 2 controlled-capacity GMM (pre-registered) | Frozen |
| 2b/2c N-scaling and EM diagnostics | Exploratory (post hoc), labelled as such |
| 2d EM attribution (pre-registered) | Completed |
| Protocol closure | Recorded (Part 6); Item 33 and Section C open |
| 3 correction experiments | Not started |

---

# Part 1 — Phase 1 (KDE)

# KDE baseline section (Phase 1 close-out) — what recursion does, and does not do, in the KDE setting

> **FROZEN (Phase 1 close-out).** No further changes to these experiments or conclusions except corrections of error. New work uses its own seed namespace and its own section.
>
> **Conclusions, with the scope in which they were tested** (replace-style recursion with reference reinjection, λ = 0.05, N = 2000, T = 59, h ∈ [0.075, 0.29], 1-D three-component target):
> 1. Population recursion is a direct reference-only scale-mixture estimator (exact for a fixed reference; with the protocol initial law, 4.6% of the T = 59 population law is residual prior P₀).
> 2. Finite-population sampling adds a major error: W₁(B,C) ≈ 0.20–0.24, matched by an iid sample of N_eff = Nλ(2−λ) = 195 at the three m tested.
> 3. Rare-mode failure is driven primarily by missing reference observations (no rare point in the retained set: 56%, 5.5% and 0% of seeds at m = 20, 100, 400). Recursion loss and spatial migration of descendants were negligible **with reinjection**.
> 4. Recursive KDE did not outperform the best reference-only estimator among those tested: a tie at m = 20 (−0.003 ± 0.017) and worse by 0.07–0.10 W₁ for m ≥ 100 (best-of-list selection on the same seeds).
>
> **Not covered:** the pure-recursion case (λ = 0, no reinjection), where rare-lineage loss is the classical mechanism, was not run for KDE, so the negative result does not extend to it. Other λ, N, T and learners are untested.


Status: results of three experiments, written as a separate section; the main protocol file is unchanged. Seed namespaces: archives 77001 (40 paired seeds, nested reference sets), 77003 (400 seeds, Experiment 2); recursion streams 88000+m, 88500+m, 99000+m, 99500+m. Setting throughout: 1-D target P0 = 0.70 N(−3,1²) + 0.27 N(2,0.8²) + 0.03 N(8,0.6²); Gaussian KDE learner; replace-style recursion; N = 2000, λ = 0.05, T = 59; with-replacement reference draw (n_ref = 100); m ∈ {20, 50, 100, 200, 400}. Scripts: `kde_common.py`, `exp_decomp_bandwidth.py`, `exp_phase1_a.py` … `exp_phase1_d.py`; outputs `*_out.txt`. Evidence tags: [A] analytic, [G] grid, [MC] Monte Carlo.

## 1. Defensible statements

1. **Population recursion is a reference-only kernel estimator [A, G].** For a fixed reference law R and h, the population operator gives q_t = R * K_t exactly, with K_t = Σ_{j<t} λ(1−λ)^j N(0,(j+1)h²) + (1−λ)^t N(0,(t+1)h²) when started from R*N_h. When started from the protocol's initial law (1−λ)P0 + λR, there is an extra term (1−λ)^{t+1} P0*N(0,(t+1)h²). At T = 59 that term has weight 0.046, so a T = 59 population law is 4.6% prior information, not purely reference. At T → ∞ the kernel is the geometric scale mixture Σ_j λ(1−λ)^j N(0,(j+1)h²), with variance h²/λ. In population the recursion therefore adds no information beyond R. It changes the smoothing kernel.
2. **Finite sampling is the dominant extra error, and it is explained by an effective sample size [MC].** The finite-N chain is a resampling-with-reinjection chain, and its lineage variance carries the factor 1/(λ(2−λ)) (Prop 3). It behaves like an iid sample of size N_eff = Nλ(2−λ) = 195 drawn from the population law.
3. **At h = 0.1 and λ = 0.05 the recursion does not destroy the rare lineage and does not move its descendants out of the rare region [MC].** Rare-region error comes from the reference sample (no rare point in the retained set) and from smoothing, not from recursion loss.
4. **Recursive C is not better than reference-only estimation except at m = 20, and there only weakly.** For m ≥ 100 it is worse by 0.07–0.10 W1.

## 2. Experiment 1: scale-mixture equivalence and finite-N discrepancy

**Equivalence [A, G].** Grid recursion A versus the direct mixture estimator B (h = 0.1, T = 59 and T = 400, both initial laws, 40 seeds per m): max W1(A,B) = 5.2e-14 and max |ΔCDF| = 1.9e-15. This is an algebraic composition identity on one grid, so it tests the closed-form weights, not the discretisation. The independent check is the analytic CDF of B (5 seeds per m): W1(A, B_exact) = 2.8e-4 to 9.1e-4, equal to W1(B_grid, B_exact), so the whole difference is binning discretisation (dx = 0.01). **Prediction W1(A,B) ≈ 0 holds to the discretisation limit of about 1e-3.**

**Finite-N discrepancy [MC]** (T = 59, h = 0.1, mean ± se over 40 paired seeds, reference-only initial law):

| m | W1(B,C) | iid floor, n = N = 2000 | iid floor, n = N_eff = 195 | W1(B,C) − floor(N_eff) |
|---|---|---|---|---|
| 20 | 0.204 ± 0.015 | 0.073 | 0.206 ± 0.016 | −0.002 ± 0.017 |
| 100 | 0.239 ± 0.018 | 0.073 | 0.206 ± 0.011 | +0.032 ± 0.018 |
| 400 | 0.238 ± 0.016 | 0.074 | 0.248 ± 0.016 | −0.010 ± 0.024 |

The discrepancy is about 3× an iid-N floor and independent of m, and is matched by the N_eff floor within error at the three m tested. With the protocol initial law the numbers are the same within error (0.183, 0.230, 0.252, 0.219, 0.233 for m = 20…400).

## 3. Experiment 2: lineage-resolved finite-N analysis (h = 0.1, 400 seeds, protocol initial law)

| m | no rare point in retained set | rare point in set, lineage lost by t = 59 | rare lineage survives | survivors: own-region retention, emp. vs pop. | finite-N excess W1(C) − W1(pop, seed-conditional) |
|---|---|---|---|---|---|
| 20 | 55.7% | 0.0% | 44.3% | 0.999 vs 0.999 (−0.0001 ± 0.0008) | +0.038 ± 0.007 |
| 100 | 5.5% | 0.25% (1 seed; ever absent 1.3%) | 94.3% | 1.000 vs 0.999 (+0.0008 ± 0.0001) | +0.078 ± 0.007 |
| 400 | 0.0% | 0.0% (ever absent 0.3%) | 100% | 0.999 vs 0.999 (+0.0005 ± 0.0003) | +0.114 ± 0.006 |

- Case (i), no rare lineage in the retained set, is the main rare-mode failure: rare-region mass error ε_R ≈ −0.95 (the mass left is the decaying remnant of the initial sample, 0.0014 at T = 59).
- Case (ii), the recursion loses a lineage it holds, is negligible at λ = 0.05 because reference reinjection reintroduces it.
- Case (iii), descendants leave the rare region, does not occur at h = 0.1: 0.1% of label-3 points lie outside R3.
- The finite-N W1 excess is similar across categories (e.g. m = 20: +0.037 for case (i), +0.039 for case (iii)). It is bulk sampling drift, not lineage loss. Because W1 is not additive, the excess grows with m only because the reference-estimation part shrinks.

## 4. Experiment 3: bandwidth neighbourhood and full trade-off (40 paired seeds, window 55–59)

C by rule, W1 (C − same-h reference-only, paired). Inflation rule h = √(ε λ v̂), v̂ from the reference sample; ε = 0.025, 0.05, 0.1, 0.2 gives h ≈ 0.10, 0.145, 0.205, 0.29.

| m | fixed 0.075 | fixed 0.1 | ε = 0.025 | ε = 0.05 | ε = 0.1 | ε = 0.2 |
|---|---|---|---|---|---|---|
| 20 | 0.700 (−0.035) | 0.693 (−0.038) | 0.699 (−0.032) | 0.703 (−0.021) | 0.732 (+0.015) | 0.803 (+0.094) |
| 50 | 0.488 (−0.002) | 0.487 (−0.001) | 0.491 (+0.003) | 0.502 (+0.018) | 0.541 (+0.062) | 0.631 (+0.158) |
| 100 | 0.413 (+0.060) | 0.412 (+0.061) | 0.415 (+0.063) | 0.428 (+0.079) | 0.467 (+0.120) | 0.552 (+0.208) |
| 200 | 0.323 (+0.072) | 0.324 (+0.074) | 0.327 (+0.076) | 0.344 (+0.094) | 0.395 (+0.147) | 0.505 (+0.255) |
| 400 | 0.270 (+0.097) | 0.273 (+0.101) | 0.275 (+0.102) | 0.295 (+0.124) | 0.353 (+0.182) | 0.462 (+0.291) |

- C is monotonically worse as ε grows; there is no interior optimum, and the best C settings are the smallest h tried. Reference-only estimators improve with larger h (best listed: ε = 0.1–0.2). The recursion-aware rule therefore does not approach the best reference-only estimator.
- Best reference-only versus best rule-based C, paired gap: m = 20: −0.003 ± 0.017 (tie); 50: +0.022 ± 0.022; 100: +0.069 ± 0.021; 200: +0.075 ± 0.019; 400: +0.099 ± 0.015. Best-of-list selection uses the same seeds on both sides and is slightly optimistic.
- Other endpoints (C, m = 100 / 200 / 400, fixed 0.1): signed ε_R = +0.100 / +0.130 / −0.085 (reference-only: +0.025 / +0.092 / +0.046); tail ℓ_τ median +0.51 / +0.72 / +0.68 and rising with h; variance ratio 1.012 / 1.041 / 0.996; label-3 absent at t = 59 in 0% of runs (m ≥ 100).
- Zero-tail frequencies: the raw reference has no tail points in 80%, 65%, 47.5%, 35% and 20% of seeds for m = 20…400. For C the fraction with model tail mass below 0.1% of the target is 0–15% (nonzero only at h = 0.075 and m ≤ 50). KDE models never have exactly zero tail mass, so ℓ_τ is reported unpadded.
- The scale-mixture reference-only estimator (population limit) at h = 0.1 gives W1 of 0.696, 0.465, 0.344 for m = 20, 50, 100, below the same-h KDE (0.731, 0.488, 0.352) and the raw reference at m = 20–50. At m = 20 the small advantage of C over reference-only KDE at the same h is therefore a kernel-shape effect.
- LOO and Silverman bandwidths are fine as reference-only estimators (m = 100: 0.355 and 0.462) but fail inside the recursion (1.04 and 1.77), because the bandwidth accumulates as h²/λ.

## 5. Limitations

- One λ (0.05), one N (2000), one horizon (T = 59), one h-family; 40 seeds in Experiments 1 and 3 and 400 in Experiment 2. Grid discretisation is about 1e-3 in W1 and ±0.003 in ε_R.
- The N_eff account matches at three m within error and is not tested at other λ or N; it is not a derivation of the W1 floor.
- Experiment 2 uses reinjection every generation; loss of the rare lineage is a no-reinjection phenomenon and is not tested here.
- Best-of-list selections in Experiment 3 are on the same seeds as the evaluation.
- Earlier figures I was given ("W1 ≈ 0.28", "reference-LOO failure") were not produced in this session and are not reconciled with these tables.


---

# Part 2 — Phase 2 (GMM) with exploratory 2b/2c

> **FROZEN (2026-11-09)** — Sections 1, 3, 4 are the Phase 2 record. Section 2 stays labelled exploratory (post hoc, 10–20 seeds). Injection wording tightened (two-thirds at m ≥ 100, not "60–70% everywhere"). Conclusion, scoped: under these conditions (1-D three-component target, N = 2000, T = 60, sklearn EM defaults, λ = 0.05 for injection), repeated finite-sample projection by GMM learners can degrade distributional fidelity; the magnitude and the rare-mode signature depend materially on the fitting procedure. No universal collapse mechanism is claimed. Follow-up: Phase 2d (`gmm_phase2d_design.md`).

# Phase 2 — controlled-capacity GMM recursion: results

Design (pre-registered before any run): `gmm_phase2_design.md`. Seed namespace 2026110800+. Scripts `gmm_phase2.py`, `gmm_phase2_analyze.py`; Phase 2b/2c are **exploratory and were added after seeing the Phase 2 results** (`gmm_phase2b.py`, `gmm_phase2c.py`). Evidence tag [MC] throughout. Phase 1 (KDE) is frozen and unchanged.

Setting: target P0 as before; N = 2000 per generation, T = 60, 60 paired seeds, m ∈ {20, 100, 400}; learners K3 (correct capacity), K2 (under-specified), BIC (K ∈ 1…6 chosen per generation); regimes R0 reference-only, R2 pure recursion (λ = 0), R3 injection (λ = 0.05), R4 fresh iid P0 data (budget control), R5 reference bootstrap. Ancestry is not tracked (GMM-generated points carry no component label).

## 1. Pre-registered results (W1, mean over seeds; late = generations 51–60)

| m | learner | R0 ref-only | R2 at t = 1 | R2 late | R3 late | R4 fresh | R2 growth (late − t=1), paired | P(rare comp) R2: t=1 → late |
|---|---|---|---|---|---|---|---|---|
| 20 | K3 | 0.709 | 0.713 | 0.880 | 0.715 | 0.129 | +0.167 ± 0.040 | 0.38 → 0.21 |
| 20 | K2 | 0.721 | 0.730 | 0.913 | 0.777 | 0.266 | +0.183 ± 0.038 | 0.05 → 0.05 |
| 20 | BIC | 0.713 | 0.715 | 0.838 | 0.750 | 0.070 | +0.123 ± 0.036 | 0.48 → 0.48 |
| 100 | K3 | 0.327 | 0.334 | 0.687 | 0.453 | 0.148 | +0.353 ± 0.033 | 0.50 → 0.11 |
| 100 | K2 | 0.393 | 0.406 | 0.744 | 0.492 | 0.268 | +0.339 ± 0.042 | 0.02 → 0.02 |
| 100 | BIC | 0.300 | 0.316 | 0.683 | 0.386 | 0.072 | +0.367 ± 0.038 | 0.83 → 0.60 |
| 400 | K3 | 0.194 | 0.248 | 0.642 | 0.394 | 0.127 | +0.393 ± 0.027 | 0.50 → 0.04 |
| 400 | K2 | 0.309 | 0.326 | 0.638 | 0.406 | 0.264 | +0.313 ± 0.028 | 0.00 → 0.00 |
| 400 | BIC | 0.156 | 0.175 | 0.521 | 0.293 | 0.072 | +0.346 ± 0.028 | 1.00 → 0.68 |

- **Exit condition.** (a) growth over generations, (b) positive slope, and (c) R2 worse than the reference bootstrap R5 and R0, each with the paired 95% CI excluding 0: **met on W1 in 9 of 9 cells; met on rare-region mass in 6 of 9 (all learners at m = 100 and 400), not at m = 20.** Slopes of W1 per generation under R2 are +0.002 to +0.003 (m = 20) and +0.005 to +0.006 (m = 100, 400), each ± 0.001; under R3 they are about +0.001 or smaller.
- **Signature.** Variance ratio stays 0.93–1.03, so it is not variance collapse. Tails and the rare component thin out: late median ℓ_τ under R2 is −1.54 (K3, m = 400) and −3.36 (BIC, m = 400); BIC rare-component presence falls 1.00 → 0.68 at m = 400 and the BIC mode count rises 3.0 → 3.7 (selected K 3.25 → 3.87), so BIC also accumulates spurious components. K2 never holds a rare component (its W1 floor is 0.26 even with fresh data).
- **Injection.** R3 removes about 63–75% of the R2 W1 growth at m ≥ 100 (and nearly all of it at m = 20); R3 − R2 late W1 = −0.09 to −0.30, CIs excluding 0 for every cell, but at m ≥ 100 R3 late W1 stays above reference-only (difference of means, R3 − R0: m = 100: +0.126, +0.099, +0.086; m = 400: +0.200, +0.097, +0.137 for K3, K2, BIC). This matches the Phase 1 KDE finding that recursion with reinjection does not beat reference-only.

## 2. Exploratory attribution (added post hoc; fewer seeds)

**N-scaling of R2** (m = 100, same references; W1 growth late − t=1, 10–20 seeds):

| learner | N = 500 | 2000 | 20,000 | 200,000 | |ε_R| growth at N = 2000 / 20,000 / 200,000 |
|---|---|---|---|---|---|
| K3 | +0.798 ± 0.145 | +0.257 ± 0.047 | +0.118 ± 0.021 | +0.067 ± 0.025 | +0.238 / +0.226 / +0.190 |
| K2 | +0.529 ± 0.092 | +0.216 ± 0.048 | +0.080 ± 0.024 | +0.038 ± 0.018 | +0.087 / +0.170 / +0.183 |
| BIC | +0.729 ± 0.139 | +0.329 ± 0.056 | +0.050 ± 0.023 | not run | +0.461 / +0.061 / – |

Growth falls with N, so most of it is finite-sample accumulation, not a fixed bias of the learner. It falls more slowly than N^(−1/2) for K3 and K2 (K3 from N = 2000 to 200,000: ×0.26 against ×0.10 predicted), and the rare-region error growth for K3 and K2 does not fall with N. K3 P(rare) reaches 0.00 at every N.

**Fitting procedure (K3, m = 100, N = 20,000, 20 seeds):**

| EM variant | W1 growth | |ε_R| growth | P(rare) t=1 → late |
|---|---|---|---|
| default (n_init = 1, k-means restart each generation) | +0.118 ± 0.021 | +0.226 ± 0.046 | 0.30 → 0.00 |
| tight (n_init = 3, tol = 1e-6) | +0.116 ± 0.022 | +0.072 ± 0.048 | 0.60 → 0.15 |
| warm-started from previous generation (tol = 1e-6) | +0.048 ± 0.026 | +0.098 ± 0.052 | 0.65 → 0.63 |

Warm starting roughly halves the W1 growth (1.8σ from zero) and retains the rare component; tighter EM alone does not. Part of the loss in Phase 2 is a re-initialisation artifact of the fitting procedure.

## 3. Statement

Under pure recursion (λ = 0), GMM learners of all three capacities show reproducible, growing degradation in W1 and, at m ≥ 100, in rare-region mass and rare-component presence. It is a finite-sample accumulation phenomenon: it shrinks with N, it is not reproduced by training on the reference alone at equal budget, and a sizeable share depends on how each generation is refit (cold start versus warm start). Reinjection at λ = 0.05 cuts the W1 growth by roughly two-thirds at m ≥ 100 without reaching reference-only fidelity.

## 4. Limitations

- Primary results at N = 2000, T = 60, λ = 0.05 for injection; one target; sklearn EM defaults (n_init = 1).
- K3 "rare presence" mixes data limits with optimisation failures: with fresh iid data (R4) only 60–72% of K3 fits place a component at the rare mode.
- The exit condition tests degradation relative to a reference bootstrap and says nothing about which mechanism causes it; Section 2 is exploratory, with 10–20 seeds, and the N = 200,000 BIC cell was not run.
- Mode count and rare presence are defined in the design; ancestry is untracked.
- No correction method was run.
- Pure KDE recursion (λ = 0) remains untested.


---

# Part 3 — Phase 2d

# Phase 2d — EM fitting-procedure comparison: results (pre-registered, `gmm_phase2d_design.md`)

K3, N = 2000, T = 60, 60 paired seeds per m, seed namespace 2026110900+. Raw tables: `gmm_phase2d_analysis.txt`. [MC]

## Pure recursion (R2), late-window W1 (generations 51–60), reference-only fit R0 = 0.353 (m = 100), 0.184 (m = 400)
| variant | m = 100 late W1 | rare presence t=1 → late | m = 400 late W1 | rare presence t=1 → late |
|---|---|---|---|---|
| V0 default | 0.688 | 0.43 → 0.02 | 0.615 | 0.55 → 0.04 |
| V1 multi-start | 0.629 | 0.55 → 0.29 | 0.581 | 0.72 → 0.26 |
| V2 warm | 0.592 | 0.53 → 0.42 | 0.540 | 0.73 → 0.49 |
| V3 hybrid | 0.635 | 0.53 → 0.44 | 0.578 | 0.75 → 0.46 |

Paired contrasts on late W1 (95% CI): V2−V0 = −0.096 [−0.188, −0.004] (m = 100), −0.075 [−0.147, −0.004] (m = 400): better, not material by the pre-set −0.05 rule on the CI. V1−V0 = −0.059 [−0.147, +0.029], −0.034 [−0.104, +0.035]: inconclusive. V3 vs V1 or V2: inconclusive in all four contrasts.

## Injection (R3, λ = 0.05)
All of V1, V2, V3 beat V0 on late W1 with CIs excluding 0 (differences −0.04 to −0.10). V3 is best at m = 400 (−0.100 [−0.135, −0.064], material; V3−V1 −0.040, V3−V2 −0.061, both CIs excluding 0). At m = 100 V1, V2, V3 are indistinguishable. R3 late minus R0: 0.05–0.15 for V1–V3 against 0.125 and 0.185 for V0; no variant fully reaches reference-only.

## Reading against the pre-registered table
- Pure recursion: V2 is better than V0, V1 is inconclusive, so the evidence leans to path continuity over local-optimum selection but does not separate them; the hybrid does not improve on either.
- No variant removes the degradation: R2 late − R0 = +0.24 to +0.43 for every variant. Fitting procedure changes the magnitude (about 15–30% of late W1 at best) and the rare-mode signature (rare presence at generation 60: 0.02–0.04 under V0, 0.26–0.49 otherwise), not the existence of degradation. The residual is not attributed to a single mechanism: finite-sample drift is one candidate; learner approximation and repeated projection are not excluded.
- Under injection the fitting procedure matters more, and combining warm start with multi-start helps at m = 400.

## Limitations
One learner (K3), one N, one T, one λ; the ±0.05 margin is judgement-based and is ratified for Phase 2d fitting-variant comparisons only (not a general W1 margin); CIs are paired-normal with 60 seeds; V0 here is the Phase 2 procedure re-run under a new namespace (its growth, +0.32 and +0.39, matches Phase 2's +0.35 and +0.39). BIC and K2 not tested.


---

# Part 4 — Mechanisms

# Mechanisms: what each phase does and does not separate

| Candidate cause of rare-mode / fidelity loss | Evidence so far | Status |
|---|---|---|
| Reference omission (the reference never contained the mode) | KDE: no rare point in the retained set in 56%, 5.5%, 0% of seeds at m = 20, 100, 400 (Phase 1). GMM: K3 on fresh iid data misses the rare component in 28–40% of seeds (R4). | Main rare-mode failure in the KDE setting [MC] |
| Finite-sample sampling | KDE: W1(B,C) ≈ 0.20–0.24, matched by an iid sample of N_eff = 195. GMM: R2 W1 growth falls with N, more slowly than N^(−1/2) for K3/K2 (Phase 2b, exploratory). | Established as a large contributor; not shown to be the only one |
| Learner approximation / optimisation path | Warm start improves pure recursion (late W1 −0.075 to −0.096, CI near 0); multi-start inconclusive; hybrid helps under injection at m = 400 (Phase 2d). | Real, partial; neither restores reference-only fidelity |
| Repeated smoothing / projection | Population KDE recursion with reinjection equals a reference-only scale-mixture estimator (no extra loss). Pure KDE recursion (λ = 0) not run. | Untested for pure recursion |

Not claimed: that the remaining GMM degradation is exclusively finite-sample drift; universal model collapse; any result outside the tested 1-D target, learners, N, and horizons. Distributional degradation is not shown to be irreversible information loss.

---

# Part 5 — Seed namespaces and files

Phase 1: 77001, 77003; recursion streams 88000+m, 88500+m, 99000+m, 99500+m. Phase 2/2b/2c: 2026110800+. Phase 2d: 2026110900+. Scripts: `kde_common.py`, `exp_phase1_*.py`, `gmm_phase2*.py`; outputs `*_out.txt`, `gmm_phase2*_analysis.txt`.

---

# Part 6 — Protocol closure

See `model_collapse_math_protocol.md` Section 16 (authoritative): stabilisation definition and floors ratified; T_fix = 60 for λ = 0 (fixed comparison horizon, not convergence); T_R = T5 stays primary (stabilisation windows are sensitivity analyses; B-specific assessment separate); Phase 3 proceeds on the existing 1-D protocol; seed plan pilot 1–100, test 1001–1300 retained; ±0.05 margin is for Phase 2d only; 2-D accumulation designs deferred.
