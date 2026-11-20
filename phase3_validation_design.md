# Phase 3 analytical validation, steps 4–5 (S-NP score; correction operator) — criteria fixed before running

Scope: validation of the implementation and characterisation of the correction operator C of Section 10 with the S-NP score. **No hypothesis (H2, AB) is tested; no correction is compared to baselines; no test seeds are used.** Seeds: validation namespace 90001+ (not pilot 1–100, not test 1001–1300). Target as in Section 2. Evidence tags: [A] analytic, [G] numerical check, [MC] Monte Carlo, [U] untested.

## V18 — S-NP score exactness [A/G]
For random reference sets (m ∈ {20, 100, 400}), s ∈ {0.1, 0.3, 1, 3, 6}, x ∈ [−15, 25]:
- (a) analytic ŝ_s(x) equals the central difference of log p_s(x) (p_s the smoothed reference KDE, via logsumexp; step 1e-5): pass if max relative error ≤ 1e-6 where |ŝ| > 1e-3, absolute ≤ 1e-8 otherwise.
- (b) denoiser D_s(x) = x + s² ŝ_s(x) equals the posterior mean Σ w_j(x) r_j: pass if max abs diff ≤ 1e-10.
- (c) DSM stationarity [MC]: J at ŝ is below J at ŝ + δ·φ (φ = smooth bounded perturbation, δ ∈ {±0.02, ±0.05}) for all four δ, same noise draws.
- (d) numerical stability: no NaN/inf for |x − r|/s up to 10³.

## V19 — s_min by reference-only LOO likelihood [G]
20 reference sets per m. Report s_min per m, whether the optimum is interior to the search grid (log-spaced 0.05–6), and the Silverman value. Flag if s_min exceeds 0.5 (then s_0 ∈ {0.5, …} cannot be below s_L, and the pilot grid of s_max must be reconsidered).

## V20 — S-NP score versus the population score of P_0 * N_s [G/MC] (descriptive)
Per m and s, mean absolute score error over regions R1, R2, R3 and the sign of the score at the rare-mode centre (x = 8) when the reference holds no rare point. No pass criterion; the result sets what the correction can possibly know about the rare mode.

## V21 — Langevin implementation check [A]
Single-Gaussian reference-law score at level s (variance v = σ_r² + s²): the discretised annealed update with step α has stationary variance v/(1 − α/(4v)). Run 10⁵ chains for 200 steps at α/v ∈ {0.05, 0.2, 0.5}: pass if the sample variance is within 2% of the formula. The step ratio of the protocol schedule is α_ℓ/s_ℓ² = ε/s_L² (constant over levels); the schedule is parametrised here by c = ε/s_L² and the stability requirement is c < 2 [A], with c ≤ 0.5 used.

## V22 — Correction operator characterisation [MC] (descriptive; conjectures C1–C3 stated in Section 10)
Setup: L = 10 geometric levels from s_0 down to s_min (s_0 ≥ s_min, so ℓ_0 ∈ {1..L}); K Langevin steps per level; c ∈ {0.05, 0.2}; s_0 ∈ {2.5, 6} ∪ {values in [s_min, 2.5]} chosen after V19; 30 reference sets (10 per m ∈ {20, 100, 400}); 4000 inputs from each of two input laws: (i) P_0 with labels, (ii) Q* = Σ_j λ(1−λ)^j P_0*N(0,(j+1)h²) with h = 0.6, λ = 0.05, labels inherited (the arm-B input in the pilot configuration).
Reported per cell: displacement |y_out − y_in| (mean, 95th percentile) by input label; label-to-region transition matrix A for input vs output and ΔA; region-3 mass and signed ε_R against P_0; tail mass; W1 of output to P_0 and to the reference-only law p_{s_min} (sample-based); mean/variance ratio.
Checks of the conjectures, with explicit criteria:
- C1: mean displacement nondecreasing in s_0 at fixed (m, K, c, input law): report the number of violations among adjacent s_0 pairs (tolerance: increase ≥ −0.01).
- C2: output depends weakly on input: |W1(out_{P0}, out_{Q*})| ≤ 0.05 W1-units and the same ordering for region-3 mass; report W1(out, p_{s_min}) next to W1(out, P_0).
- C3: for s_0 ≤ 1, off-diagonal |ΔA| (label k ending in a region other than its input region) ≤ 0.005 for every label.
Outcomes are recorded as met / not met; neither choice of regime is selected here (that is the pilot's job, using seeds 1–100).

## Not covered
Steps 1–3, 6–7 of the Phase 3 validation list; S-NN; baselines; full recursion with correction; power; seeds 1–100 and 1001–1300.
