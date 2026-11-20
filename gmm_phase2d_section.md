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
