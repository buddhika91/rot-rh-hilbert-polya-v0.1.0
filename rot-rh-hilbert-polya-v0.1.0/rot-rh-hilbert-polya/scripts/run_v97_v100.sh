#!/usr/bin/env bash
set -euo pipefail
REPO="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
RESEARCH="$REPO/research/legacy_pipeline"
WORK="$REPO/runs/v97_v100"
mkdir -p "$WORK"
export PYTHONPATH="$RESEARCH${PYTHONPATH:+:$PYTHONPATH}"
cd "$WORK"
python "$RESEARCH/rot_rh_full_rot_hp_pipeline_v97.py" \
  --base-module rot_rh_pure_mangoldt_prime_phase_dirac_jacobi_operator_audit \
  --L-list 2000,4000,8000,16000 --depths 24,32,48 --tau 0.02 \
  --tail-factor 18 --scan-points 96 --moment-order 6 --identity-max 50 \
  --identity-step 0.25 --xi-dps 100 --zero-scan-max 120 --zero-scan-step 0.10 \
  --zero-match-count 32 --maximum-tail-phase-bound 2e-7 --maximum-root-residual 1e-8 \
  --prefix rot_rh_full_rot_hp_pipeline_v97
python "$RESEARCH/rot_rh_weyl_tail_completion_v98.py" \
  --v97-prefix rot_rh_full_rot_hp_pipeline_v97 --depths 24,32,48 --L-final 16000 \
  --tail-moment-order 16 --tail-order-check 12 --xi-dps 100 --prefix rot_rh_weyl_tail_completion_v98
python "$RESEARCH/rot_rh_tail_midpoint_s1_v99.py" \
  --v97-prefix rot_rh_full_rot_hp_pipeline_v97 --base-module rot_rh_pure_mangoldt_prime_phase_dirac_jacobi_operator_audit \
  --L 16000 --depth 48 --tail-factor 18 --xi-dps 60 --prefix rot_rh_tail_midpoint_s1_v99
python "$RESEARCH/rot_rh_full_arithmetic_tail_v100_1.py" \
  --v97-prefix rot_rh_full_rot_hp_pipeline_v97 --base-module rot_rh_pure_mangoldt_prime_phase_dirac_jacobi_operator_audit \
  --L 16000 --depths 32,48 --tail-factor 18 --tail-moment-order 10 --post-moment-order 4 \
  --identity-max 50 --identity-step 0.5 --xi-dps 60 --prefix rot_rh_full_arithmetic_tail_v100_1
