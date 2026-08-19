$ErrorActionPreference = "Stop"
$Repo = Split-Path -Parent $PSScriptRoot
$Research = Join-Path $Repo "research\legacy_pipeline"
$Work = Join-Path $Repo "runs\v102_1"
New-Item -ItemType Directory -Force -Path $Work | Out-Null
$env:PYTHONPATH = $Research
Push-Location $Work
try {
  python (Join-Path $Research "rot_rh_pnt_finite_part_full_pipeline_v102_1.py") `
    --base-module rot_rh_pure_mangoldt_prime_phase_dirac_jacobi_operator_audit `
    --v97-prefix (Join-Path $Repo "runs\v97_v100\rot_rh_full_rot_hp_pipeline_v97") `
    --L-list 2000,4000,8000,16000 --depths 24,32,48 `
    --tau 0.02 --tail-factor 18 --laguerre-q 384 --laguerre-qcheck 768 `
    --scan-points 96 --tail-moment-order 10 --post-moment-order 4 `
    --identity-max 50 --identity-step 0.5 --xi-dps 60 `
    --prefix rot_rh_pnt_finite_part_full_pipeline_v102_1
} finally { Pop-Location }
