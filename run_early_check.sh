#!/bin/bash
cd /home/claude
for pair in "exp_phase1_a.py exp_phase1_a_out.txt" "exp_phase1_b.py exp_phase1_b_out.txt" "exp_phase1_d.py exp_phase1_d_out.txt" "exp_phase1_c.py exp_phase1_c_out.txt" "validate_rev4.py validate_rev4_out.txt" "validate_rev4b.py validate_rev4b_out.txt" "validate_rev5_horizon.py validate_rev5_horizon_out.txt" "validate_rev5_sim.py validate_rev5_sim_out.txt" "validate_rev5_sim2.py validate_rev5_sim2_out.txt" "validate_rev5_traj.py validate_rev5_traj_out.txt" "phase3_val_a.py phase3_val_a_out.txt" "phase3_val_a2.py phase3_val_a2_out.txt" "phase3_step2_lam0.py phase3_step2_out.txt" "exp_decomp_bandwidth.py exp_decomp_bandwidth_out.txt"; do
 set -- $pair; s=$(date +%s); timeout 900 python3 $1 > /tmp/early_$1.txt 2>&1; rc=$?; e=$(( $(date +%s)-s ))
 if diff -q /tmp/early_$1.txt $2 >/dev/null; then echo "$1: IDENTICAL (rc=$rc, ${e}s)"; else echo "$1: DIFFERS rc=$rc ${e}s, $(diff /tmp/early_$1.txt $2 | grep -c '^[<>]') differing lines"; fi
done
echo ALLDONE
