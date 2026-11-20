# Reproducing the model-collapse programme (Phase 1 → Section 30)

## Environment (as run)
Python 3.13.18, numpy 2.5.3, scipy 1.18.1, scikit-learn 1.9.1, mpmath 1.4.1 (see `requirements_used.txt`). Linux x86-64, 2 vCPU (all multiprocessing uses `Pool(2)`; BLAS threads forced to 1 inside the scripts). No GPU. No network at run time. No git history: protocol versions are kept as `rev5_pre_*_backup.md`.
Run everything from this directory (scripts read/write the current directory).

## Verification levels and what was tested (2026-11-10)
| Level | Result |
|---|---|
| Same environment, one seed re-run from scratch, compared with saved pickle | **max abs difference 0.0**: stage-1 pilot seed 1; Sections 25 (seed 201), 26 (31), 27 (61), 28 (301), 30 (1001); GMM Phase 2d (m=100, seed 0) |
| Different environment (numpy 2.1.3, scipy 1.14.1, scikit-learn 1.6.1) | Sections 25–30 KDE blocks: **0.0**. GMM Phase 2d: V0, V1, V2 exact; **V3 (multi-start EM) differs** |
| Earlier phases re-run in full against saved outputs | Phase 1 `exp_phase1_a/b/c/d`, `exp_decomp_bandwidth`, `validate_rev4`, `rev4b`, `rev5_horizon`, `rev5_sim`, `rev5_sim2`, `rev5_traj`, `phase3_val_a`, `phase3_step2_lam0`: identical (numbers; `validate_rev5_sim` differs only in printed wall-clock times). `phase3_val_a2` prints the output of `phase3_val_a` first because `phase3_val_a` has no `__main__` guard (its own part is identical). Not re-run in full: `phase3_val_b.py` (analysis re-run on its pickle: identical), GMM Phase 2/2b/2c (analysis re-run on pickles: identical), full stage-1 and Phase 2d (one seed/job each). |
| Analysis scripts re-run on saved pickles | identical text output for Sections 25, 26, 27, 28, 30, stage 1, CDF diagnostic, Phase 2/2b/2c/2d, V21–V22 |
| Numbers in protocol result sections vs analysis outputs | automatic check (`audit_numbers.py`): all non-matching decimals were explained (derived quantities, constants, rounding) |

## Run order (all Phase 3 blocks; wall times on 2 cores)
`./reproduce.sh full` runs, in order: stage-1 pilot (`phase3_pilot.py 1 30 stage1`, ~70 min), CDF diagnostic (~9 min), Section 25 pilot (~5 min), Block 1 (seeds 31–60, ~2 min), Block 2 (61–90, ~20 s), Block 3 (301–330, ~1 s), confirmatory H1 (1001–1300, ~1 min), then analyses and `export_csv.py`. `./reproduce.sh verify` (~3 min) re-runs one seed per block and all analyses. Earlier phases: `exp_phase1_*.py`, `exp_decomp_bandwidth.py`, `validate_rev*.py`, `phase3_val_*.py`, `phase3_step2_lam0.py`; GMM `gmm_phase2.py`, `gmm_phase2b.py`, `gmm_phase2c.py`, `gmm_phase2d.py` (7, 4, 3, 46 min) each followed by its `*_analyze.py`.

## Seed registry
Pilot stage 1 / CDF diagnostic 1–30; Block 1 31–60; Block 2 61–90; reserve 91–100; Section 25 pilot 201–230; reserved for deferred mixing diagnostic 231–260; Block 3 301–330; test 1001–1300 used for H1 only (Section 30). Every seed gets independent `SeedSequence([seed, tag, …])` streams; seed ranges are disjoint, so no stream is shared across blocks. Phase 1/2 use fixed integer seeds in the scripts (77001–99500; `BASE2 = 2026110900`).

## Known issues
1. **Library-version sensitivity of EM multi-start (Phase 2d V3).** Under scikit-learn 1.6.1 the V3 trajectories differ in 12/12 checked seeds (m = 100: R3 mean W1 at T 0.392 saved vs 0.326; R2 0.679 vs 0.637; paired |diff| 0.07–0.10). V0–V2 are unaffected. V3-related contrasts in Phase 2d should be read as version-specific; pin scikit-learn 1.9.1 to reproduce them.
2. CLI arguments are undocumented in the scripts (`phase3_pilot.py a b tag`, `phase3_pilot_analyze.py tag`, `phase3_s2[678].py a b file`); see `reproduce.sh`.
3. Results are stored as Python pickles (version- and trust-sensitive); `export_s2*_per_seed.csv` and `export_s30_per_seed.csv` give per-seed results as plain CSV for Sections 25–30 only. Phase 1–2 and stage-1 pickles have no CSV export.
4. Analysis bootstrap CIs depend on seed ordering; the analysis scripts sort by seed, so outputs are stable.
5. Wall-clock times and compute counters in `phase3_step6_out.txt` vary by machine.
6. `pilot_stage1_audit_2seeds.pkl` is the audit trail for the tail-guard amendment (Section 20.B).
7. Timing-dependent behaviour: runs are resumable (`phase3_pilot.py`, `phase3_s25_pilot.py`, s26–s28); they were interrupted and resumed by sandbox restarts, and resumed results equal uninterrupted ones because every seed is a self-contained job.
