## 21. Pilot stage 1 results (seeds 1–30, m = 100; screening only) — `phase3_pilot.py`, `phase3_pilot_analyze.py`, `phase3_pilot_stage1_analysis.txt`, `phase3_pilot_stage1.pkl`
Rules applied exactly as frozen in Sections 20.A and 20.B (tail ceiling 0.0248). Seeds 1–30 only; no configuration beyond stage 1 has been run; the test seeds are untouched. Evidence [MC]. A first launch was lost to a sandbox restart and rerun from scratch with identical code (resumable); the two-seed audit file `pilot_stage1_audit_2seeds.pkl` is from the rerun's first checkpoint.

**Outcome under the frozen rule: all 12 configurations are eliminated** (lower 95% Wilson bound of at least one single-guardrail failure rate exceeds 0.10 for every configuration). No stage 2 has been started.

| Guardrail | Failure-rate range over the 12 configurations | Reference-only KDE (R_h) on the same seeds | Uncorrected C |
|---|---|---|---|
| Variance (relative to C) | 0 for 11; 0.17 for s₀ = 2.5, c = 0.05, 10×10 | n/a | n/a |
| Tail ≤ 0.0248 | 0.36–1.00 (best: s₀ = 1.0, c = 0.2, 20×20 at 0.36 [0.20, 0.55]) | 0.03 | 1.00 |
| Region 3 in [0.015, 0.060] | 0.13–0.67 | 0.33 | 0.33 |
| Original tail ≤ 0.009 (audit) | 0.96–1.00 | 0.73 | 1.00 |

**What the failures are.**
- *Region-3 guard:* reference-only KDE fails it on 33% of seeds, and B's ρ₃ tracks R_h's (correlation 0.89–0.98 across seeds, 2/30 references contain no rare point). The guard therefore measures how many rare points a 100-point reference happened to hold, not the stability of the correction. A 10% failure rule cannot be met by any estimator that uses these references, including reference-only KDE.
- *Tail guard:* informative, not an artefact of the reference. R_h passes (tail 0.0129 mean, 3% failures), C fails every seed (mean 0.119), and B sits in between: mean tail 0.0205–0.0616 across configurations, about 1.75–2.0× the reference-only KDE tail at the better schedules, with the smoothing floor at 0.0124 and the unsmoothed target at 0.004533. At best B's tail is 36% failures at the 0.0248 ceiling.

**Endpoint observations (screening only; seeds 1–30; not confirmatory).**
- B reduces W₁ by 0.43–0.61 relative to C at the schedules with s₀ ≤ 1 or c = 0.2 (se ≈ 0.03–0.04); the exceptions are s₀ = 2.5, c = 0.05, 10×10 (−0.09).
- No configuration has a lower W₁ than the reference-only KDE at h = 0.6: B − R_h = +0.068 to +0.565 (R_h W₁ = 0.338, LOO-bandwidth reference-only 0.348). The best B configurations (20×20, s₀ = 1.0 or 2.5, c = 0.2) are within about +0.07 of R_h.
- Variance ratio at the R window is 1.09–1.30 for the better configurations (C: 1.80).
- Rare-lineage own-region retention A₃₃ for B is 0.11–0.73 against 0.85–0.86 for C (higher with s₀ = s_min: 0.59–0.73; 0.11–0.32 at s₀ = 1.0 or 2.5, with c = 0.2). Non-rare inflow into region 3: 0.010–0.066 for B against 0.026 for C; J₃ (lineage-based net flow): +0.002 to +0.042 for B against +0.022 for C. (Averaged over seeds where the quantity is defined; reference-set SD not shown.)
- Compute: 10×10 about 10 s per run, 20×20 about 40 s (100 and 400 updates per correction).
- s_min across the 30 references: median 0.53 (0.27–1.09). The s₀ lower-bound rule raised the nominal 1.0 in 5 seeds, and 2.5 in none; those five seeds have fewer than 3 distinct s₀ values (5 de-duplicated settings).

**Reading.** In this schedule family the correction removes most of the uncorrected recursion's error but does not beat reference-only KDE on W₁ and leaves the tail about twice as heavy as reference-only KDE. These are stage-1 observations on 30 seeds and say nothing yet about H2, which requires the confirmatory design and all baselines.
