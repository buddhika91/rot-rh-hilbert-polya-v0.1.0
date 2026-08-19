import json
from pathlib import Path

def test_reported_release_status_is_not_proof():
    root=Path(__file__).resolve().parents[1]
    d=json.loads((root/'results/reported_metrics.json').read_text())
    assert d['v100_1']['rh_proof'] is False
    assert d['v100_1']['xi_zeta_known_zeros_pre_freeze'] is False
    assert d['v98']['weyl_completed_det_symmetric_rms'] < d['v98']['raw_det_symmetric_rms']
