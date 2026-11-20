# Phase 3 analytical validation — steps 4–5 (S-NP score, correction operator)

Criteria: `phase3_validation_design.md` (fixed before running). Scripts `phase3_val_a.py`, `phase3_val_a2.py` (diagnostic), `phase3_val_b.py`, `phase3_val_b_analyze.py`; outputs `phase3_val_a_out.txt`, `phase3_val_a2_out.txt`, `phase3_val_b_analysis.txt`. Validation seed namespace 90001+ (pilot 1–100 and test 1001–1300 untouched). No hypothesis tested.

## Results
| Check | Outcome |
|---|---|
| V18(a) analytic score vs finite difference | **Criterion not met as written** (max rel. error 2.0e-5 vs ≤ 1e-6). Post hoc diagnostic: against a 40-digit mpmath evaluation the analytic score agrees to 4.4e-16 (s = 0.1); a 4th-order stencil gives 5e-10 for s ≥ 1 and 2e-4 only at s = 0.1. The miss is finite-difference error at small s, not the formula. [A/G] |
| V18(b) D_s = posterior mean | Met: max abs diff 5.8e-12 |
| V18(c) DSM stationarity | Met: J rises for all four perturbations (0.5730 → 0.5732–0.5740) [MC] |
| V18(d) extreme separation | Met: finite at \|x − r\|/s = 10⁵ |
| V19 s_min (LOO) | Interior optimum in 60/60 sets. Median s_min 0.79 (m = 20; range 0.25–2.46), 0.51 (m = 100), 0.35 (m = 400); all well below Silverman (1.57, 1.23, 0.92). Flag raised: s_min exceeds the planned lower s_0 in some m = 20 sets, so s_0 is taken as max(s_0, s_min) |
| V20 S-NP vs population score of P0*N_s | Descriptive. At s = 1 the mean score error in the rare region is 2.8 (m = 20), 1.3 (m = 100), 0.12 (m = 400). When the reference has no rare point (5/10, 2/10, 0/10 refs) the score at x = 8 is strongly negative (−5.3 and −5.2 at s = 1; −0.67 at s = 3, against a true −0.38 at s = 3): the score pulls samples away from the missing mode |
| V21 Langevin stationary variance | Met: sample variance within 0.71% of v/(1 − α/4v) at α/v = 0.05, 0.2, 0.5 [A/MC] |
| V22 C1 (displacement increases with s_0) | **Met**: 0 violations in 36 adjacent pairs |
| V22 C2 (output depends weakly on input, \|ΔW1\| ≤ 0.05) | **Partly met**: met for large s_0 (2.5, 6.0) at c = 0.2 (7 of 9 such cells; m = 100 and 400 at s_0 = 2.5 pass at 0.038 and 0.045), not met at s_0 ≤ 1 or at c = 0.05 (up to 0.56) |
| V22 C3 (no inter-mode transport at s_0 ≤ 1, max off-diagonal \|ΔA\| ≤ 0.005) | **Not met** in 24 of 24 cells (range 0.006–0.134). The closest cell is m = 400, s_0 = s_min, c = 0.05 (0.0059) |

## What the characterisation shows (n = 2000 inputs, 10 references per m; L = 10 levels, K = 10 steps per level)
- Fed a perfect P0 sample (input W1 to P0 = 0.08–0.09), the correction returns something at about reference quality: output W1 to P0 is 0.58–0.92 (m = 20), 0.19–0.59 (m = 100) and 0.11–0.60 (m = 400) depending on (s_0, c); at s_0 = s_min, c = 0.2 it is 0.69, 0.22, 0.13. It cannot return more than the reference holds.
- Fed the arm-B input Q* (input W1 to P0 = 0.95–0.99), it reduces W1 to 0.37–0.77 at best (s_0 ∈ {s_min, 1.0}, c = 0.2), but the output still depends on the input: at m = 400 the Q* output W1 is 0.44 against 0.13 for the same operator on a P0 input. In the pilot configuration the correction therefore recovers part of the gap, not all of it; no comparison with reference-only or with baselines has been made.
- c = 0.05 with K = 10, L = 10 is under-annealed at large s_0: at s_0 = 6 the tail mass is 0.12–0.17 (P0: 0.0045), region-3 mass 0.09–0.12 (P0: 0.030), variance 16–21 (P0: 8.4). c = 0.2 gives 0.02, 0.05–0.07 and 11–12. Choosing (L, K, c) is therefore a first-order pilot decision, not a detail.
- Inter-mode transport is present at every setting tested, increasing with s_0 and c, and decreasing with m (e.g. summed off-diagonal |ΔA| at s_0 = s_min, c = 0.05: 0.17 at m = 20, 0.04 at m = 100, 0.011 at m = 400). The correction should not be assumed to be within-mode reshaping at any s_0 ≥ s_min.

## Limits
Single-step operator on fixed input laws; the full recursion with correction was not run. One schedule family (L = 10, K = 10, geometric levels, c ∈ {0.05, 0.2}); s_0 grid relative to s_min. 10 references per m; 2000 inputs; off-diagonal ΔA has Monte Carlo noise of order 0.003, which is relevant to the 0.005 C3 threshold. C3's threshold was set by me and is not a ratified margin. S-NN, baselines, power and steps 1–3, 6–7 of the validation list are not done.
