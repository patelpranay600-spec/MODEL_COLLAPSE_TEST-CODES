# Protocol closure (Phase 1–2d) — decisions recorded 2026-11-09

## A. Ratified decisions
1. **Stabilisation definition (ratified).** For endpoint E with reference level E* (the uncorrected-operator limit) and anchor E0, t_stab is the first generation t such that for every generation s ≥ t up to the horizon, |E_s − E*| ≤ max(0.05·|E* − E0|, ε_abs(E)), **and** this holds over a window of 10 consecutive generations (sustained crossing, all-subsequent condition retained).
2. **Numerical floors ε_abs (adopted as the proposed implementation).**

| Endpoint | ε_abs |
|---|---|
| W1 | 0.001 |
| Rare-region mass | 0.0001 |
| Variance ratio | 0.001 |
| Log-tail measure | 0.01 |

3. **λ = 0 arms.** T5 is undefined (no contraction). A separate **fixed horizon** is used: T_fix = 60 generations, matching Phase 2 and the Phase 2d design; the primary comparison for λ = 0 is generation-matched at T_fix. *Proposed by the assistant; needs your confirmation.* Stabilisation is not reported for λ = 0 unless the trajectory meets the rule inside T_fix.
4. **Equivalence margin.** ±0.05 on late-window W1 is ratified for **Phase 2d fitting-variant comparisons only**. It is not a general W1 margin; Phase 3 margins are set separately before the test seeds.

## B. Not decided (stays open)
- **Item 33, R-horizon rule** T_R = max(T5, t_stab,e_R, t_stab,e_τ, t_stab,W1). The ratification covers the stabilisation definition and floors, not this rule. Until you ratify it, the primary R window stays at T5 with the horizon-limited flags of Section 7.3.
- The R-window width (5 generations, [T_R−4, T_R]) and the 10-generation stability window are different quantities (the first is the evaluation window, the second is the stabilisation criterion).
- Item 31 (P-H2c as a fraction of the C-to-oracle-R gap) and the lineage-loss equivalence feasibility.

## C. Discrepancy to resolve before Phase 3 step 2
The Phase 3 roadmap refers to an "earlier frozen design" with a **2-D Gaussian KDE**, rare-mode probability 3% (1% alternative via pilot), **accumulate-subsample (primary) and full accumulation (secondary)**, a 2 × 2 factorial, an **analytic reference-trained score**, EM-GMM K = 2 as secondary comparison, and **100 pilot + 300 fresh test seeds**. The protocol on disk (Revision 4 + addendum) records a **1-D** three-component target (rare weight 0.03), replace-style recursion D_t = (1−λ)·synthetic + λ·reference with a fixed reference budget, and explicitly defers the 2-D extension "until the 1-D mechanism is understood" (Section 15, item 19). It contains no accumulate-subsample/full-accumulation arms and no 100/300 seed plan. I have no record of that design being frozen. Phase 3 analytical validation therefore needs either (a) your confirmation that the 2-D accumulate-subsample design replaces the 1-D replace-style recursion for Phase 3, with its spec supplied, or (b) Phase 3 run on the existing 1-D replace-style protocol. The conditions A–D in the roadmap (noise-assisted, correction, uncorrected, noise + correction) map onto the existing factorial regardless.

## D. Evidence classes for conclusions
- Phase 1 (KDE, replace-style with reinjection, 1-D): [A]/[G] for the scale-mixture identity and lineage variance; [MC] for finite-N and bandwidth results.
- Phase 2 (GMM, pre-registered): [MC], 60 paired seeds.
- Phase 2b/2c: [MC], **exploratory, post hoc**, 10–20 seeds.
- Phase 2d: [MC], pre-registered, 60 paired seeds.
