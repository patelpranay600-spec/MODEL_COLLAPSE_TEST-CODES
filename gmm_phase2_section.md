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
