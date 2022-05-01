import numpy as np


def calc_fadums_from_m_and_n(mval, nval):
    """
    Computes the influence factor for stress at the corner of a foundation according to Fadum's Chart.

    Equations from Bowles (1996)

    Parameters
    ----------
    mval
    nval

    Returns
    -------

    """
    i = np.where(mval ** 2 * nval ** 2 > mval ** 2 + nval ** 2 + 1, (2 * mval * nval * (mval ** 2 + nval ** 2 + 2) / (
            np.sqrt(mval ** 2 + nval ** 2 + 1) * (mval ** 2 + nval ** 2 + 1 + mval ** 2 * nval ** 2)) + np.arctan(
        2 * mval * nval * np.sqrt(mval ** 2 + nval ** 2 + 1) / (mval ** 2 + nval ** 2 + 1 - mval ** 2 * nval ** 2)) + np.pi) / (
                         4 * np.pi), (2 * mval * nval * (mval ** 2 + nval ** 2 + 2) / (
            np.sqrt(mval ** 2 + nval ** 2 + 1) * (mval ** 2 + nval ** 2 + 1 + mval ** 2 * nval ** 2)) + np.arctan(
        2 * mval * nval * np.sqrt(mval ** 2 + nval ** 2 + 1) / (mval ** 2 + nval ** 2 + 1 - mval ** 2 * nval ** 2))) / (4 * np.pi))
    return i


def calc_stress_under_corner(length, width, depth, stress):
    m = length / depth
    n = width / depth
    i = calc_fadums_from_m_and_n(m, n)
    return stress * i


def calc_stress_under_centre(length, width, depth, stress):
    m = length / depth / 2
    n = width / depth / 2
    i = calc_fadums_from_m_and_n(m, n) * 4
    return stress * i
