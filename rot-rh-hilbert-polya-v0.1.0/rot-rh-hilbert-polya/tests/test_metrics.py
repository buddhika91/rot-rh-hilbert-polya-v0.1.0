from rot_rh_hp.metrics import relative_l2, symmetric_rms

def test_metrics_zero_on_identity():
    a=[1.0,2.0,3.0]
    assert relative_l2(a,a) == 0.0
    assert symmetric_rms(a,a) == 0.0
