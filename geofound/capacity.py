"""
Created on Jan 17, 2014

@author: maximmillen
"""

import numpy as np

from geofound.output import log
from geofound.exceptions import DesignError
from geofound import models


def vesics_capacity(sl, fd, h_l=0, h_b=0, vertical_load=1, slope=0, base_tilt=0, verbose=0):
    """
    Calculates the foundation capacity according Vesics(1975)
    #Gunaratne, Manjriker. 2006. "Spread Footings: Analysis and Design."
    Ref: http://geo.cv.nctu.edu.tw/foundation/download/
                            BearingCapacityOfFoundations.pdf
    :param sl: Soil object
    :param fd: Foundation object
    :param h_l: Horizontal load parallel to length
    :param h_b: Horizontal load parallel to width
    :param vertical_load: Vertical load
    :param slope: ground slope
    :param base_tilt: The slope of the underside of the foundation
    :param verbose: verbosity
    :return: ultimate bearing stress
    """
    area_foundation = fd.length * fd.width
    c_a = 0.6 - 1.0 * sl.cohesion

    horizontal_load = np.sqrt(h_l ** 2 + h_b ** 2)

    fd.nq_factor = ((np.tan(np.pi / 4 + sl.phi_r / 2)) ** 2 * np.exp(np.pi * np.tan(sl.phi_r)))
    if sl.phi_r == 0:
        fd.nc_factor = 5.14
    else:
        fd.nc_factor = (fd.nq_factor - 1) / np.tan(sl.phi_r)
    fd.ng_factor = 2.0 * (fd.nq_factor + 1) * np.tan(sl.phi_r)

    # shape factors:
    s_c = 1.0 + fd.nq_factor / fd.nc_factor * fd.width / fd.length
    s_q = 1 + fd.width / fd.length * np.tan(sl.phi_r)
    s_g = max(1.0 - 0.4 * fd.width / fd.length, 0.6)  # add limit of 0.6 based on Vesic
    # depth factors:
    if fd.depth / fd.width > 1:
        k = np.arctan(fd.depth / fd.width)
    else:
        k = fd.depth / fd.width
    d_c = 1 + 0.4 * k
    d_q = 1 + 2 * np.tan(sl.phi_r) * (1 - np.sin(sl.phi_r)) ** 2 * k
    d_g = 1.0

    # load inclination factors
    m__b = (2.0 + fd.width / fd.length) / (1 + fd.width / fd.length)
    m_l = (2.0 + fd.length / fd.width) / (1 + fd.length / fd.width)
    m = np.sqrt(m__b ** 2 + m_l ** 2)

    if sl.phi_r == 0:
        i_q = 1.0
        i_g = 1.0
    else:
        i_q = (1.0 - horizontal_load / (vertical_load + area_foundation *
                                        c_a / np.tan(sl.phi_r))) ** m
        i_g = (1.0 - horizontal_load / (vertical_load + area_foundation *
                                        c_a / np.tan(sl.phi_r))) ** (m + 1)
    i_c = i_q - (1 - i_q) / (fd.nq_factor - 1)
    check_i_c = 1 - m * horizontal_load / (area_foundation * c_a * fd.nc_factor)
    if abs(check_i_c - i_c) / i_c > 0.001:
        raise DesignError

    # ground slope factors:
    if sl.phi_r == 0:
        # g_c = slope / 5.14
        g_c = i_q
    else:
        g_c = i_q - (1 - i_q) / (5.14 * np.tan(sl.phi_r))
    g_q = (1.0 - np.tan(slope)) ** 2
    g_g = g_q

    # tilted base factors
    if sl.phi_r == 0:
        b_c = g_c
    else:
        b_c = 1 - 2 * base_tilt / (5.14 * np.tan(sl.phi_r))
    b_q = (1.0 - base_tilt * np.tan(sl.phi_r)) ** 2
    b_g = b_q

    # stress at footing base:
    q_d = sl.unit_dry_weight * fd.depth

    if verbose:
        log("Nc: ", fd.nc_factor)
        log("N_qV: ", fd.nq_factor)
        log("Ng: ", fd.ng_factor)
        log("s_c: ", s_c)
        log("s_q: ", s_q)
        log("s_g: ", s_g)
        log("d_c: ", d_c)
        log("d_q: ", d_q)
        log("d_g: ", d_g)
        log("i_c: ", i_c)
        log("i_q: ", i_q)
        log("i_g: ", i_g)
        log("g_c: ", g_c)
        log("g_q: ", g_q)
        log("g_g: ", g_g)
        log("b_c: ", b_c)
        log("b_q: ", b_q)
        log("b_g: ", b_g)
        log("q_d: ", q_d)

    # Capacity
    fd.q_ult = (sl.cohesion * fd.nc_factor * s_c * d_c * i_c * g_c * b_c +
                q_d * fd.nq_factor * s_q * d_q * i_q * g_q * b_q +
                0.5 * fd.width * sl.unit_dry_weight *
                fd.ng_factor * s_g * d_g * i_g * g_g * b_g)

    if verbose:
        log("qult: ", fd.q_ult)
    return fd.q_ult


def terzaghi_capacity(sl, fd, round_footing=False, verbose=0):
    """
    Calculates the foundation capacity according Terzaghi (1943)
    Ref: http://geo.cv.nctu.edu.tw/foundation/
    download/BearingCapacityOfFoundations.pdf
    :param sl: Soil object
    :param fd: Foundation object
    :param round_footing: if True, then foundation is round
    :param verbose: verbosity
    :return: ultimate bearing stress
    Note: the shape factor of 1.3 is used for aspect ratio > 6
    """

    a02 = ((np.exp(np.pi * (0.75 - sl.phi / 360) * np.tan(sl.phi_r))) ** 2)
    a0_check = (np.exp((270 - sl.phi) / 180 * np.pi * np.tan(sl.phi_r)))
    if (a02 - a0_check) / a02 > 0.001:
        raise DesignError

    fd.nq_factor = (a02 / (2 * (np.cos((45 + sl.phi / 2) * np.pi / 180)) ** 2))
    fd.ng_factor = (2 * (fd.nq_factor + 1) * np.tan(sl.phi_r) / (1 + 0.4 * np.sin(4 * sl.phi_r)))
    if sl.phi_r == 0:
        fd.nc_factor = 5.7
    else:
        fd.nc_factor = (fd.nq_factor - 1) / np.tan(sl.phi_r)

    # shape factors:
    if round_footing:
        s_c = 1.3
        s_g = 0.6
    elif fd.length / fd.width < 5:
        s_c = 1.3
        s_g = 0.8
    else:
        s_c = 1.0
        s_g = 1.0
    s_q = 1.0

    # stress at footing base:
    q_d = sl.unit_dry_weight * fd.depth

    # Capacity
    fd.q_ult = (sl.cohesion * fd.nc_factor * s_c + q_d * fd.nq_factor * s_q + 0.5 * fd.width *
                sl.unit_dry_weight * fd.ng_factor * s_g)

    if verbose:
        log("Nc: ", fd.nc_factor)
        log("Nq: ", fd.nq_factor)
        log("Ng: ", fd.ng_factor)
        log("s_c: ", s_c)
        log("s_q: ", s_q)
        log("s_g: ", s_g)
        log("qult: ", fd.q_ult)
    return fd.q_ult


def hansen_capacity(sl, fd, h_l=0, h_b=0, vertical_load=1, slope=0, base_tilt=0, verbose=0):
    """
    Calculates the foundation capacity according Hansen (1970)
    Ref: http://bestengineeringprojects.com/civil-projects/
    hansens-bearing-capacity-theory/

    :param sl: Soil object
    :param fd: Foundation object
    :param h_l: Horizontal load parallel to length
    :param h_b: Horizontal load parallel to width
    :param vertical_load: Vertical load
    :param slope: ground slope
    :param base_tilt: The slope of the underside of the foundation
    :param verbose: verbosity
    :return: ultimate bearing stress
    """

    area_foundation = fd.length * fd.width
    horizontal_load = np.sqrt(h_l ** 2 + h_b ** 2)
    c_a = 0.6 - 1.0 * sl.cohesion

    fd.nq_factor = ((np.tan(np.pi / 4 + sl.phi_r / 2)) ** 2 * np.exp(np.pi * np.tan(sl.phi_r)))
    if sl.phi_r == 0:
        fd.nc_factor = 5.14
    else:
        fd.nc_factor = (fd.nq_factor - 1) / np.tan(sl.phi_r)
    fd.ng_factor = 1.5 * (fd.nq_factor - 1) * np.tan(sl.phi_r)

    # shape factors
    if sl.phi_r == 0:
        s_c = 0.2 * fd.width / fd.length
    else:
        s_c = 1.0 + fd.nq_factor / fd.nc_factor * fd.width / fd.length

    s_q = 1.0 + fd.width / fd.length * np.sin(sl.phi_r)
    s_g = 1.0 - 0.4 * fd.width / fd.length

    # depth factors:
    if fd.depth / fd.width > 1:
        k = np.arctan(fd.depth / fd.width)
    else:
        k = fd.depth / fd.width
    d_c = 1 + 0.4 * k
    if sl.phi == 0:
        d_c = 0.4 * k
    d_q = 1 + 2 * np.tan(sl.phi_r) * (1 - np.sin(sl.phi_r)) ** 2 * k
    d_g = 1.0

    # incline load factors:
    if sl.phi_r == 0:
        i_q = 1.0
        i_c = 0.5 - 0.5 * np.sqrt(1 - horizontal_load / area_foundation * c_a)
        i_g = 1.0
    else:
        i_q = ((1.0 - 0.5 * horizontal_load /
                (vertical_load + area_foundation * c_a / np.tan(sl.phi_r))) ** 5)
        i_c = i_q - (1 - i_q) / (fd.nq_factor - 1)
        i_g = ((1 - (0.7 * horizontal_load) /
                (vertical_load + area_foundation * c_a / np.tan(sl.phi_r))) ** 5)

    # slope factors:
    if sl.phi_r == 0:
        g_c = (slope / np.pi * 180) / 147
    else:
        g_c = 1.0 - (slope / np.pi * 180) / 147
    g_q = 1 - 0.5 * np.tan(slope) ** 5
    g_g = g_q

    # base tilt factors:
    if sl.phi_r == 0:
        b_c = (base_tilt / np.pi * 180) / 147
    else:
        b_c = 1.0 - (base_tilt / np.pi * 180) / 147
    b_q = (np.exp(-0.0349 * (base_tilt / np.pi * 180) * np.tan(sl.phi_r)))
    b_g = (np.exp(-0.0471 * (base_tilt / np.pi * 180) * np.tan(sl.phi_r)))

    if verbose:
        log("Nc: ", fd.nc_factor)
        log("Nq: ", fd.nq_factor)
        log("Ng: ", fd.ng_factor)
        log("s_c: ", s_c)
        log("s_q: ", s_q)
        log("s_g: ", s_g)
        log("d_c: ", d_c)
        log("d_q: ", d_q)
        log("d_g: ", d_g)
        log("i_c: ", i_c)
        log("i_q: ", i_q)
        log("i_g: ", i_g)
        log("g_c: ", g_c)
        log("g_q: ", g_q)
        log("g_g: ", g_g)
        log("b_c: ", b_c)
        log("b_q: ", b_q)
        log("b_g: ", b_g)

    # stress at footing base:
    q_d = sl.unit_dry_weight * fd.depth

    # Capacity
    if sl.phi_r == 0:
        fd.q_ult = (sl.cohesion * fd.nc_factor *
                    (1 + s_c + d_c - i_c - g_c - b_c) + q_d)
    else:
        fd.q_ult = (sl.cohesion * fd.nc_factor *
                    s_c * d_c * i_c * g_c * b_c +
                    q_d * fd.nq_factor * s_q * d_q * i_q * g_q * b_q +
                    0.5 * fd.width * sl.unit_dry_weight *
                    fd.ng_factor * s_g * d_g * i_g * g_g * b_g)


def meyerhoff_capacity(sl, fd, h_l=0, h_b=0, vertical_load=1, verbose=0):
    """
    Calculates the foundation capacity according Meyerhoff (1963)
    http://www.engs-comp.com/meyerhof/index.shtml

    :param sl: Soil object
    :param fd: Foundation object
    :param h_l: Horizontal load parallel to length
    :param h_b: Horizontal load parallel to width
    :param vertical_load: Vertical load
    :param verbose: verbosity
    :return: ultimate bearing stress
    """

    horizontal_load = np.sqrt(h_l ** 2 + h_b ** 2)

    fd.nq_factor = ((np.tan(np.pi / 4 + sl.phi_r / 2)) ** 2 *
                 np.exp(np.pi * np.tan(sl.phi_r)))
    if sl.phi_r == 0:
        fd.nc_factor = 5.14
    else:
        fd.nc_factor = (fd.nq_factor - 1) / np.tan(sl.phi_r)
    fd.ng_factor = (fd.nq_factor - 1) * np.tan(1.4 * sl.phi_r)

    if verbose:
        log("Nc: ", fd.nc_factor)
        log("Nq: ", fd.nq_factor)
        log("Ng: ", fd.ng_factor)

    kp = (np.tan(np.pi / 4 + sl.phi_r / 2)) ** 2
    # shape factors
    s_c = 1 + 0.2 * kp * fd.width / fd.length
    if sl.phi > 10:
        s_q = 1.0 + 0.1 * kp * fd.width / fd.length
    else:
        s_q = 1.0
    s_g = s_q

    # depth factors
    d_c = 1 + 0.2 * np.sqrt(kp) * fd.depth / fd.width
    if sl.phi > 10:
        d_q = 1 + 0.1 * np.sqrt(kp) * fd.depth / fd.width
    else:
        d_q = 1.0
    d_g = d_q

    # inclination factors:
    theta_load = np.arctan(horizontal_load / vertical_load)
    i_c = (1 - theta_load / (np.pi * 0.5)) ** 2
    i_q = i_c
    if sl.phi > 0:
        i_g = (1 - theta_load / sl.phi_r) ** 2
    else:
        i_g = 0

    # stress at footing base:
    q_d = sl.unit_dry_weight * fd.depth

    # Capacity
    fd.q_ult = (sl.cohesion * fd.nc_factor * s_c * d_c * i_c +
                q_d * fd.nq_factor * s_q * d_q * i_q +
                0.5 * fd.width * sl.unit_dry_weight *
                fd.ng_factor * s_g * d_g * i_g)
    return fd.q_ult


def nzs_vm4_capacity(sl, fd, h_l=0, h_b=0, vertical_load=1, slope=0, verbose=0, **kwargs):
    """
    calculates the capacity according to
     Appendix B verification method 4 of the NZ building code

    :param sl: Soil object
    :param fd: Foundation object
    :param h_l: Horizontal load parallel to length
    :param h_b: Horizontal load parallel to width
    :param vertical_load: Vertical load
    :param slope: ground slope
    :param verbose: verbosity
    :return: ultimate bearing stress
    """
    # Need to make adjustments if sand  has DR<40% or
    # clay has liquidity indices greater than 0.7

    horizontal_load = np.sqrt(h_l ** 2 + h_b ** 2)

    h_eff_b = kwargs.get("h_eff_b", 0)
    h_eff_l = kwargs.get("h_eff_l", 0)
    loc_v_l = kwargs.get("loc_v_l", fd.length / 2)
    loc_v_b = kwargs.get("loc_v_b", fd.width / 2)

    ecc_b = h_b * h_eff_b / vertical_load
    ecc_l = h_l * h_eff_l / vertical_load

    width_eff = min(fd.width, 2 *
                    (loc_v_b + ecc_b), 2 * (fd.width - loc_v_b - ecc_b))
    length_eff = min(fd.length, 2 *
                     (loc_v_l + ecc_l), 2 * (fd.length - loc_v_l - ecc_l))
    area_foundation = length_eff * width_eff

    # check para 3.4.1
    if width_eff / 2 < fd.width / 6:
        raise DesignError("failed on eccentricity")

    # LOAD FACTORS:
    fd.nq_factor = ((np.tan(np.pi / 4 + sl.phi_r / 2)) ** 2 * np.exp(np.pi * np.tan(sl.phi_r)))
    if sl.phi_r == 0:
        fd.nc_factor = 5.14
    else:
        fd.nc_factor = (fd.nq_factor - 1) / np.tan(sl.phi_r)
    fd.ng_factor = 2.0 * (fd.nq_factor - 1) * np.tan(sl.phi_r)

    # shape factors:
    s_c = 1.0 + fd.nq_factor / fd.nc_factor * width_eff / length_eff
    s_q = 1 + width_eff / length_eff * np.tan(sl.phi_r)
    s_g = max(1.0 - 0.4 * width_eff / length_eff, 0.6)  # add limit of 0.6 based on Vesics

    # depth factors:
    if fd.depth / width_eff > 1:
        k = np.arctan(fd.depth / width_eff)
    else:
        k = fd.depth / width_eff
    if sl.phi_r == 0:
        d_c = 1 + 0.4 * k
        d_q = 1.0
    else:
        d_q = (1 + 2 * np.tan(sl.phi_r) *
               (1 - np.sin(sl.phi_r)) ** 2 * k)
        d_c = d_q - (1 - d_q) / (fd.nq_factor * np.tan(sl.phi_r))
    d_g = 1.0

    # load inclination factors:
    if sl.phi_r == 0:
        i_c = 0.5 * (1 + np.sqrt(1 - horizontal_load / (area_foundation * sl.cohesion)))
        i_q = 1.0
        i_g = 1.0
    else:
        if h_b == 0:
            i_q = 1 - horizontal_load / (vertical_load + area_foundation * sl.cohesion /
                                         np.tan(sl.phi_r))
            i_g = i_q
        elif h_b > 0 and h_l == 0:
            i_q = ((1 - 0.7 * horizontal_load / (vertical_load + area_foundation * sl.cohesion /
                                                 np.tan(sl.phi_r))) ** 3)
            i_g = ((1 - horizontal_load / (vertical_load + area_foundation * sl.cohesion /
                                           np.tan(sl.phi_r))) ** 3)
        else:
            raise DesignError("not setup for bi-directional loading")
        i_c = (i_q * fd.nq_factor - 1) / (fd.nq_factor - 1)

    # ground slope factors:
    g_c = 1 - slope * (1.0 - fd.depth / (2 * width_eff)) / 150
    g_q = (1 - np.tan(slope * (1 - fd.depth / (2 * width_eff)))) ** 2
    g_g = g_q

    # stress at footing base:
    q_d = sl.unit_dry_weight * fd.depth

    if verbose:
        log("Nc: ", fd.nc_factor)
        log("Nq: ", fd.nq_factor)
        log("Ng: ", fd.ng_factor)
        log("H: ", horizontal_load)
        log("s_c: ", s_c)
        log("s_q: ", s_q)
        log("s_g: ", s_g)
        log("d_c: ", d_c)
        log("d_q: ", d_q)
        log("d_g: ", d_g)
        log("i_c: ", i_c)
        log("i_q: ", i_q)
        log("i_g: ", i_g)
        log("g_c: ", g_c)
        log("g_q: ", g_q)
        log("g_g: ", g_g)

    # Capacity
    fd.q_ult = (sl.cohesion * fd.nc_factor * s_c * d_c * i_c * g_c +
                q_d * fd.nq_factor * s_q * d_q * i_q * g_q +
                0.5 * width_eff * sl.unit_dry_weight *
                fd.ng_factor * s_g * d_g * i_g * g_g)
    if verbose:
        log("q_ult: ", fd.q_ult)
    return fd.q_ult


def salgado_capacity(sl, fd, h_l=0, h_b=0, vertical_load=1, verbose=0, **kwargs):
    """
    calculates the capacity according to
     THe Engineering of Foundations textbook by Salgado

    :param sl: Soil object
    :param fd: Foundation object
    :param h_l: Horizontal load parallel to length
    :param h_b: Horizontal load parallel to width
    :param vertical_load: Vertical load
    :param verbose: verbosity
    :return: ultimate bearing stress
    """
    # Need to make adjustments if sand  has DR<40% or
    # clay has liquidity indices greater than 0.7

    h_eff_b = kwargs.get("h_eff_b", 0)
    h_eff_l = kwargs.get("h_eff_l", 0)
    loc_v_l = kwargs.get("loc_v_l", fd.length / 2)
    loc_v_b = kwargs.get("loc_v_b", fd.width / 2)

    ecc_b = h_b * h_eff_b / vertical_load
    ecc_l = h_l * h_eff_l / vertical_load

    width_eff = min(fd.width, 2 * (loc_v_b + ecc_b), 2 * (fd.width - loc_v_b - ecc_b))
    length_eff = min(fd.length, 2 * (loc_v_l + ecc_l), 2 * (fd.length - loc_v_l - ecc_l))

    # check para 3.4.1
    if width_eff / 2 < fd.width / 6:
        DesignError("failed on eccentricity")

    # LOAD FACTORS:
    fd.nq_factor = np.exp(np.pi * np.tan(sl.phi_r)) * (1 + np.sin(sl.phi_r)) / (1 - np.sin(sl.phi_r))
    fd.ng_factor = 1.5 * (fd.nq_factor - 1) * np.tan(sl.phi_r)
    # fd.ng_factor = (fd.nq_factor - 1) * np.tan(1.32 * sl.phi_r)
    if sl.phi_r == 0:
        fd.nc_factor = 5.14
    else:
        fd.nc_factor = (fd.nq_factor - 1) / np.tan(sl.phi_r)

    # shape factors:
    s_q = 1 + (width_eff / length_eff) * np.tan(sl.phi_r)
    s_g = max(1 - 0.4 * width_eff / length_eff, 0.6)
    s_c = 1.0

    # depth factors:
    d_q = 1 + 2 * np.tan(sl.phi_r) * (1 - np.sin(sl.phi_r)) ** 2 * fd.depth / width_eff
    d_g = 1.0
    d_c = 1.0

    # stress at footing base:
    q_d = sl.unit_dry_weight * fd.depth

    if verbose:
        log("width_eff: ", width_eff)
        log("length_eff: ", length_eff)
        log("Nc: ", fd.nc_factor)
        log("Nq: ", fd.nq_factor)
        log("Ng: ", fd.ng_factor)
        log("s_c: ", s_c)
        log("s_q: ", s_q)
        log("s_g: ", s_g)
        log("d_c: ", d_c)
        log("d_q: ", d_q)
        log("d_g: ", d_g)
        log("q_d: ", q_d)

    # Capacity
    fd.q_ult = (sl.cohesion * fd.nc_factor * s_c * d_c +
                q_d * fd.nq_factor * s_q * d_q +
                0.5 * width_eff * sl.unit_dry_weight *
                fd.ng_factor * s_g * d_g)

    if verbose:
        log("qult: ", fd.q_ult)
    return fd.q_ult


def size_footing(sl, vertical_load, fos=1.0, length_to_width=1.0, verbose=0, **kwargs):
    """
    Determine the size of a footing given an aspect ratio and a load
    :param sl: Soil object
    :param vertical_load: The applied load to the foundation
    :param fos: The target factor of safety
    :param length_to_width: The desired length to width ratio of the foundation
    :param verbose: verbosity
    :return: a Foundation object
    """
    method = kwargs.get("method", 'vesics')
    depth_to_width = kwargs.get("depth_to_width", 0)
    depth = kwargs.get("depth", 0)
    use_depth_to_width = 0
    if not depth:
        use_depth_to_width = 1

    # Find approximate size
    fd = models.Foundation()
    fd.width = .5  # start with B=1.0m
    for i in range(50):
        fd.length = length_to_width * fd.width
        if use_depth_to_width:
            fd.depth = depth_to_width * fd.width
        method_selector(sl, fd, method)
        q = fd.q_ult

        bearing_capacity = q * fd.length * fd.width
        fs_actual = bearing_capacity / vertical_load

        if fs_actual < fos:
            # Need to increase foundation sizes
            fd.width += 0.5
        else:
            if verbose:
                log("fs_actual: ", fs_actual)
                log("fd.width: ", fd.width)
            break

    # at this stage the current size should be too big
    width_array = []
    fs_array = []
    for j in range(11):
        width_array.append(fd.width)
        fd.length = length_to_width * fd.width
        if use_depth_to_width:
            fd.depth = depth_to_width * fd.width
        method_selector(sl, fd, method)
        q = fd.q_ult

        capacity = q * fd.length * fd.width

        fs_array.append(capacity / vertical_load)

        fd.width = fd.width - 0.5 / 10

    # search the array until FS satisfied:
    if verbose:
        log("reqFS: ", fos)
        log("width array: \n", width_array)
        log("FS array: \n", fs_array)

    for fs in range(len(fs_array)):
        if fs_array[fs] < fos:
            fd.width = width_array[fs - 1]
            fd.length = length_to_width * fd.width
            if use_depth_to_width:
                fd.depth = depth_to_width * fd.width
            method_selector(sl, fd, method)
            break
        if fs == len(fs_array) - 1:
            DesignError("No suitable foundation sizes could be determined!")

    return fd


def method_selector(sl, fd, method, **kwargs):
    """
    Calculates the bearing capacity of a foundation on soil using the specified method.
    :param sl: Soil Object
    :param fd: Foundation Object
    :param method: Method
    :param kwargs:
    :return:
    """

    if method == 'vesics':
        vesics_capacity(sl, fd, **kwargs)
    elif method == 'nzs':
        nzs_vm4_capacity(sl, fd, **kwargs)
    elif method == 'terzaghi':
        terzaghi_capacity(sl, fd, **kwargs)
    elif method == 'hansen':
        hansen_capacity(sl, fd, **kwargs)
    elif method == 'meyerhoff':
        meyerhoff_capacity(sl, fd, **kwargs)
    elif method == 'salgado':
        salgado_capacity(sl, fd, **kwargs)

available_methods = {
    "vesics": vesics_capacity,
    "nzs": nzs_vm4_capacity,
    "terzaghi": terzaghi_capacity,
    "hansen": hansen_capacity,
    "meyerhoff": meyerhoff_capacity,
    "salgado": salgado_capacity
}


def capacity_from_spt():
    pass
    # Capacity:
    # #Determine the undrained strength (Su) from CPT (Bowles, 1995):
    # Su = (qc - sigma_p_v0)/N_k
    # where N_k ranges from 15-30 depending on
    # consolidation and plastic index.
    # N_k = 13+5.5/50*PI
    # #Determine allowable bearing capacity from SPT (Parry, 1977):
    # q_n_all = 30N_55*(s/25.4) #or check section 4.3.1
    # #or convert to phi
    # phi = 25 + 28 * (N_55 / q)**0.5
