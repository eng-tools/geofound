import numpy as np

import geofound as gf


def test_strip_stress():
    depth = 4.0  # m
    q_b = 150.0
    expected_q = 46.0
    stress = gf.fadums_chart.calc_stress_under_strip_footing(width=2, depth=depth, stress=q_b)
    assert np.isclose(stress, expected_q, rtol=0.03), (stress, expected_q)

