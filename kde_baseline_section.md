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
