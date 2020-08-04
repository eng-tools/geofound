import numpy as np


def interp2d(x, xf, f):
    """
    Can interpolate a table to get an array of values in 2D

    Parameters
    ----------
    x: array_like
        1d array of values to be interpolated
    xf: 1d array of values
    f: array_like
        2d array of function values size=(len(xf), n)

    Returns
    -------

    Examples
    --------
    >>> f = np.array([[0, 0, 0],
    >>>              [0, 1, 4],
    >>>              [2, 6, 2],
    >>>              [10, 10, 10]
    >>>              ])
    >>> xf = np.array([0, 1, 2, 3])

    >>> x = np.array([0.5, 1, 2.2, 2.5])
    >>> f_interp = interp2d(x, xf, f)
    >>> print(f_interp[0][0])
    0.0
    >>> print(f_interp[0][1])
    0.5
    >>> print(f_interp[0][2])
    2.0
    """
    ind = np.argmin(np.abs(x[:, np.newaxis] - xf), axis=1)
    x_ind = xf[ind]
    ind0 = np.where(x_ind > x, ind - 1, ind)
    ind1 = np.where(x_ind > x, ind, ind + 1)
    ind0 = np.clip(ind0, 0, None)
    ind1 = np.clip(ind1, None, len(xf) - 1)
    f0 = f[ind0]
    f1 = f[ind1]
    a0 = xf[ind0]
    a1 = xf[ind1]
    denom = (a1 - a0)
    denom_adj = np.clip(denom, 1e-10, None)  # to avoid divide by zero warning
    s0 = np.where(denom > 0, (x - a0) / denom_adj, 1)  # if denom less than 0, then out of bounds
    s1 = 1 - s0
    return s1[:, np.newaxis] * f0 + s0[:, np.newaxis] * f1


def interp_twoways(x0, y0, xs, ys, fs):
    f_interp = interp2d(np.array([y0]), ys, fs)[0]
    return np.interp(x0, xs, f_interp)


def get_kz_gazetas_v_lte_0p4(a0, lob):
    lobs = np.array([1.2, 4, 6, 10, 1000])
    a0s = np.linspace(0, 2, 10)
    kzs = np.array([
        [1, 1, 0.9885, 0.9677, 0.9216, 0.8732, 0.8147, 0.7504, 0.6791, 0.6001],
        [1, 1.012, 1.032, 1.042, 1.023, 0.986, 0.9348, 0.8665, 0.7861, 0.6985],
        [1, 1.069, 1.142, 1.199, 1.196, 1.135, 1.043, 0.9441, 0.8508, 0.7575],
        [1, 1.124, 1.225, 1.275, 1.252, 1.184, 1.094, 0.9916, 0.8906, 0.8016],
        [1, 1.231, 1.348, 1.356, 1.309, 1.233, 1.143, 1.047, 0.9523, 0.8579],
    ])
    return interp_twoways(a0, lob, a0s, lobs, kzs)


def get_kz_gazetas_v_gt_0p4(a0, lob):
    lobs = np.array([1, 6, 1000])
    a0s = np.linspace(0, 2, 10)
    kzs = np.array([
        [1, 1.007, 0.9804, 0.9235, 0.853, 0.7575, 0.6508, 0.5295, 0.3859, 0.2301],
        [1, 1.03, 1.043, 0.9886, 0.8232, 0.567, 0.2662, -0.1174, -0.5485, -0.5485],
        [1, 1.227, 1.266, 1.117, 0.805, 0.3359, -0.2942, -0.7138, -0.7138, -0.7138],
        ])
    return interp_twoways(a0, lob, a0s, lobs, kzs)


def get_ky_gazetas(a0, lob):
    lobs = np.array([1, 2, 4, 6, 10, 1000])
    a0s = np.linspace(0, 2, 10)
    kys = np.array([
        [1, 1.033, 1.035, 1.014, 0.9868, 0.9596, 0.9346, 0.9174, 0.9025, 0.8892],
        [1, 1.05, 1.087, 1.09, 1.09, 1.088, 1.083, 1.081, 1.083, 1.092],
        [1, 1.099, 1.165, 1.214, 1.248, 1.272, 1.28, 1.275, 1.264, 1.263],
        [1, 1.142, 1.245, 1.317, 1.366, 1.4, 1.418, 1.423, 1.423, 1.421],
        [1, 1.193, 1.344, 1.47, 1.539, 1.564, 1.572, 1.577, 1.575, 1.573],
        [1, 1.237, 1.485, 1.631, 1.699, 1.723, 1.726, 1.729, 1.728, 1.725],
        ])
    return interp_twoways(a0, lob, a0s, lobs, kys)


def get_cz_gazetas_v_lte_0p4(a0, lob):
    lobs = np.array([1, 2, 4, 6, 10, 1000])
    a0s = np.linspace(0, 2, 10)
    czs = np.array([
        [0.9178, 0.9252, 0.9327, 0.9402, 0.9478, 0.9555, 0.9633, 0.9711, 0.9856, 1.001],
        [0.9716, 0.9784, 0.9852, 0.9921, 0.9988, 1.005, 1.011, 1.017, 1.025, 1.037],
        [1.109, 1.102, 1.095, 1.088, 1.08, 1.07, 1.061, 1.057, 1.059, 1.061],
        [1.486, 1.29, 1.174, 1.114, 1.085, 1.072, 1.063, 1.06, 1.06, 1.061],
        [1.773, 1.416, 1.252, 1.176, 1.143, 1.121, 1.106, 1.096, 1.092, 1.088],
        [2.003, 1.579, 1.344, 1.235, 1.183, 1.152, 1.132, 1.122, 1.117, 1.115],
    ])
    return interp_twoways(a0, lob, a0s, lobs, czs)


def get_czf_gazetas_v_gt_0p4(a0, lob):
    lobs = np.array([1, 2, 4, 6, 1000])
    a0s = np.linspace(0, 2, 10)
    czfs = np.array([
        [1.038, 1.03, 1.022, 1.015, 1.009, 1.006, 1.003, 1.003, 1.004, 1.004],
        [1.038, 1.03, 1.022, 1.015, 1.009, 1.006, 1.003, 1.003, 1.004, 1.004],
        [1.077, 1.061, 1.047, 1.035, 1.028, 1.023, 1.019, 1.017, 1.017, 1.017],
        [1.12, 1.096, 1.077, 1.062, 1.052, 1.045, 1.04, 1.036, 1.034, 1.036],
        [1.223, 1.158, 1.11, 1.078, 1.062, 1.053, 1.047, 1.046, 1.046, 1.049],
        ])
    return interp_twoways(a0, lob, a0s, lobs, czfs)


if __name__ == '__main__':
    print(get_kz_gazetas_v_lt_0p4(1, 1.5))
