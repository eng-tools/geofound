import numpy as np
from geofound import tables_of_dyn_stiffness_coefficients as tdsc

def test_interp2d():
    x = np.linspace(1, 10, 3)
    xf = np.linspace(0, 22, 5)
    f = np.arange(len(xf))[:, np.newaxis] * np.ones((len(xf), 10))
    f_interp = tdsc.interp2d(x, xf, f)
    assert np.isclose(f_interp[0][0], (x[0] - xf[0]) / (xf[1] - xf[0])), (f_interp[0][0], (x[0] - xf[0]) / (xf[1] - xf[0]))
    assert np.isclose(f_interp[1][0], 1), (f_interp[0][0], 1)
    assert len(f_interp) == 3
    assert len(f_interp[0]) == 10


def test_interp2d_2():
    f = np.array([[0, 0, 0],  # 0
                  [0, 1, 4],  # 5
                  [2, 6, 2],  # 10
                  [10, 10, 10]  # 30
                  ])
    xf = np.array([0, 1, 2, 3])

    x = np.array([0.5, 1, 2.2, 2.5])
    f_interp = tdsc.interp2d(x, xf, f)
    print(f_interp)
    assert f_interp[0][0] == 0
    assert f_interp[0][1] == 0.5
    assert f_interp[0][2] == 2.0
    assert f_interp[1][0] == f[1][0]
    assert f_interp[1][1] == f[1][1]
    assert f_interp[1][2] == f[1][2]
    assert np.isclose(f_interp[2][0], f[2][0] + 8 * 0.2)
    assert np.isclose(f_interp[3][2], 6.)


def test_interp2d_at_edge():
    f = np.array([[0, 0, 0],  # 0
                  [10, 10, 10]  # 30
                  ])
    xf = np.array([0, 3])

    x = np.array([0.0, 3.0])
    f_interp = tdsc.interp2d(x, xf, f)
    print(f_interp)
    assert f_interp[0][0] == 0
    assert f_interp[1][0] == 10.


def test_get_kz_gazetas_v_lt_0p4():
    assert np.isclose(tdsc.get_kz_gazetas_v_lt_0p4(1, 1.5), 0.93756)