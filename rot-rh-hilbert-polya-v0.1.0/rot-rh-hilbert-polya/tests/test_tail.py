import numpy as np
from rot_rh_hp.tail import weyl_tail_moments, midpoint_moments, finite_det, completed_det


def test_weyl_tail_is_positive_and_decays_for_relevant_cutoff():
    w = weyl_tail_moments(140.0, 5)
    assert np.all(w > 0)
    assert np.all(w[1:] < w[:-1])


def test_midpoint_is_negative():
    m = midpoint_moments(140.0, 5)
    assert np.all(m < 0)


def test_finite_det_vanishes_at_a_root():
    roots = np.array([2.0, 3.0, 5.0])
    d = finite_det(roots, [0.0, 2.0])
    assert d[0] == 1.0
    assert abs(d[1]) < 1e-15


def test_zero_tail_leaves_determinant_unchanged():
    roots = [2.0, 3.0]
    t = np.array([0.0, 0.25, 0.5])
    raw = finite_det(roots, t)
    got = completed_det(raw, t, np.zeros(4))
    assert np.allclose(raw, got)
