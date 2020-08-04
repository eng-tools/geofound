import numpy as np

from geofound.output import log
from geofound.exceptions import DesignError
from geofound import models
import sfsimodels as sm


def capacity_vesic_1975(sl, fd, h_l=0, h_b=0, vertical_load=1, slope=0, base_tilt=0, verbose=0, gwl=1e6, **kwargs):
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

    if not kwargs.get("disable_requires", False):
        models.check_required(sl, ["phi_r", "cohesion", "unit_dry_weight"])
        models.check_required(fd, ["length", "width", "depth"])
    if fd.length > fd.width:  # TODO: deal with plane strain
        fd_length = fd.length
        fd_width = fd.width
    else:
        fd_length = fd.width
        fd_width = fd.length
    area_foundation = fd_length * fd_width
    c_a = 0.6 - 1.0 * sl.cohesion

    horizontal_load = np.sqrt(h_l ** 2 + h_b ** 2)

    fd.nq_factor = ((np.tan(np.pi / 4 + sl.phi_r / 2)) ** 2 * np.exp(np.pi * np.tan(sl.phi_r)))
    if sl.phi_r == 0:
        fd.nc_factor = 5.14
    else:
        fd.nc_factor = (fd.nq_factor - 1) / np.tan(sl.phi_r)
    fd.ng_factor = 2.0 * (fd.nq_factor + 1) * np.tan(sl.phi_r)

    # shape factors:
    s_c = 1.0 + fd.nq_factor / fd.nc_factor * fd_width / fd_length
    s_q = 1 + fd_width / fd_length * np.tan(sl.phi_r)
    s_g = max(1.0 - 0.4 * fd_width / fd_length, 0.6)  # add limit of 0.6 based on Vesic
    # depth factors:
    if fd.depth / fd_width > 1:
        k = np.arctan(fd.depth / fd_width)
    else:
        k = fd.depth / fd_width
    d_c = 1 + 0.4 * k
    d_q = 1 + 2 * np.tan(sl.phi_r) * (1 - np.sin(sl.phi_r)) ** 2 * k
    d_g = 1.0

    # load inclination factors
    m__b = (2.0 + fd_width / fd_length) / (1 + fd_width / fd_length)
    m_l = (2.0 + fd_length / fd_width) / (1 + fd_length / fd_width)
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

    if gwl == 0:
        q_d = sl.unit_eff_weight * fd.depth
        unit_weight = sl.unit_bouy_weight

    elif gwl > 0 and gwl < fd.depth:
        q_d = (sl.unit_dry_weight * gwl) + (sl.unit_bouy_weight * (fd.depth - gwl))
        unit_weight = sl.unit_bouy_weight

    elif gwl >= fd.depth and gwl <= fd.depth + fd_width:
        sl.average_unit_bouy_weight = sl.unit_bouy_weight + (
            ((gwl - fd.depth) / fd_width) * (sl.unit_dry_weight - sl.unit_bouy_weight))
        q_d = sl.unit_dry_weight * fd.depth
        unit_weight = sl.average_unit_bouy_weight

    elif gwl > fd.depth + fd_width:
        q_d = sl.unit_dry_weight * fd.depth
        unit_weight = sl.unit_dry_weight

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
                0.5 * fd_width * unit_weight *
                fd.ng_factor * s_g * d_g * i_g * g_g * b_g)

    if verbose:
        log("qult: ", fd.q_ult)
    return fd.q_ult

def capacity_vesics_1975(sl, fd, h_l=0, h_b=0, vertical_load=1, slope=0, base_tilt=0, verbose=0, gwl=1e6, **kwargs):
    return capacity_vesic_1975(sl, fd, h_l=0, h_b=0, vertical_load=1, slope=0, base_tilt=0, verbose=0, gwl=1e6, **kwargs)


def capacity_terzaghi_1943(sl, fd, round_footing=False, verbose=0, **kwargs):
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
    if not kwargs.get("disable_requires", False):
        models.check_required(sl, ["phi_r", "cohesion", "unit_dry_weight"])
        models.check_required(fd, ["length", "width", "depth"])
    if fd.length > fd.width:  # TODO: deal with plane strain
        fd_length = fd.length
        fd_width = fd.width
    else:
        fd_length = fd.width
        fd_width = fd.length

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
    elif fd_length / fd_width < 5:
        s_c = 1.3
        s_g = 0.8
    else:
        s_c = 1.0
        s_g = 1.0
    s_q = 1.0

    # stress at footing base:
    q_d = sl.unit_dry_weight * fd.depth

    # Capacity
    fd.q_ult = (sl.cohesion * fd.nc_factor * s_c + q_d * fd.nq_factor * s_q + 0.5 * fd_width *
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


def capacity_hansen_1970(sl, fd, h_l=0, h_b=0, vertical_load=1, slope=0, base_tilt=0, verbose=0, **kwargs):
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
    if not kwargs.get("disable_requires", False):
        models.check_required(sl, ["phi_r", "cohesion", "unit_dry_weight"])
        models.check_required(fd, ["length", "width", "depth"])

    if fd.length > fd.width:  # TODO: deal with plane strain
        fd_length = fd.length
        fd_width = fd.width
    else:
        fd_length = fd.width
        fd_width = fd.length
    area_foundation = fd_length * fd_width
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
        s_c = 0.2 * fd_width / fd_length
    else:
        s_c = 1.0 + fd.nq_factor / fd.nc_factor * fd_width / fd_length

    s_q = 1.0 + fd_width / fd_length * np.sin(sl.phi_r)
    s_g = 1.0 - 0.4 * fd_width / fd_length

    # depth factors:
    if fd.depth / fd_width > 1:
        k = np.arctan(fd.depth / fd_width)
    else:
        k = fd.depth / fd_width
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
                    0.5 * fd_width * sl.unit_dry_weight *
                    fd.ng_factor * s_g * d_g * i_g * g_g * b_g)


def capacity_meyerhof_1963(sl, fd, gwl=1e6, h_l=0, h_b=0, vertical_load=1, verbose=0, **kwargs):
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
    if not kwargs.get("disable_requires", False):
        models.check_required(sl, ["phi_r", "cohesion", "unit_dry_weight"])
        models.check_required(fd, ["length", "width", "depth"])
    horizontal_load = np.sqrt(h_l ** 2 + h_b ** 2)

    if fd.length > fd.width:  # TODO: deal with plane strain
        fd_length = fd.length
        fd_width = fd.width
    else:
        fd_length = fd.width
        fd_width = fd.length

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
    s_c = 1 + 0.2 * kp * fd_width / fd_length
    if sl.phi > 10:
        s_q = 1.0 + 0.1 * kp * fd_width / fd_length
    else:
        s_q = 1.0
    s_g = s_q

    # depth factors
    d_c = 1 + 0.2 * np.sqrt(kp) * fd.depth / fd_width
    if sl.phi > 10:
        d_q = 1 + 0.1 * np.sqrt(kp) * fd.depth / fd_width
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

    if gwl == 0:
        q_d = sl.unit_bouy_weight * fd.depth
        unit_weight = sl.unit_bouy_weight

    elif gwl > 0 and gwl < fd.depth:
        q_d = (sl.unit_dry_weight * gwl) + (sl.unit_bouy_weight * (fd.depth - gwl))
        unit_weight = sl.unit_bouy_weight

    elif gwl >= fd.depth and gwl <= fd.depth + fd_width:
        sl.average_unit_bouy_weight = sl.unit_bouy_weight + (
        ((gwl - fd.depth) / fd_width) * (sl.unit_dry_weight - sl.unit_bouy_weight))
        q_d = sl.unit_dry_weight * fd.depth
        unit_weight = sl.average_unit_bouy_weight

    elif gwl > fd.depth + fd_width:
        q_d = sl.unit_dry_weight * fd.depth
        unit_weight = sl.unit_dry_weight

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
        log("q_d: ", q_d)

    # Capacity
    fd.q_ult = (sl.cohesion * fd.nc_factor * s_c * d_c * i_c +
                q_d * fd.nq_factor * s_q * d_q * i_q +
                0.5 * fd_width * unit_weight *
                fd.ng_factor * s_g * d_g * i_g)
    return fd.q_ult


def capacity_nzs_vm4_2011(sl, fd, h_l=0, h_b=0, vertical_load=1, slope=0, verbose=0, **kwargs):
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

    if not kwargs.get("disable_requires", False):
        models.check_required(sl, ["phi_r", "cohesion", "unit_dry_weight"])
        models.check_required(fd, ["length", "width", "depth"])

    horizontal_load = np.sqrt(h_l ** 2 + h_b ** 2)

    v_l = kwargs.get("loc_v_l", fd.length / 2)  # given in fd coordinates
    v_b = kwargs.get("loc_v_b", fd.width / 2)
    h_eff_b = kwargs.get("h_eff_b", 0)
    h_eff_l = kwargs.get("h_eff_l", 0)

    if fd.length > fd.width:  # TODO: deal with plane strain
        fd_length = fd.length
        fd_width = fd.width
        loc_v_l = v_l
        loc_v_b = v_b
        loc_h_b = h_b
        loc_h_l = h_l
    else:
        fd_length = fd.width
        fd_width = fd.length
        loc_v_l = v_b
        loc_v_b = v_l
        loc_h_b = h_l
        loc_h_l = h_b

    ecc_b = loc_h_b * h_eff_b / vertical_load
    ecc_l = loc_h_l * h_eff_l / vertical_load

    width_eff = min(fd_width, 2 *
                    (loc_v_b + ecc_b), 2 * (fd_width - loc_v_b - ecc_b))
    length_eff = min(fd_length, 2 *
                     (loc_v_l + ecc_l), 2 * (fd_length - loc_v_l - ecc_l))
    area_foundation = length_eff * width_eff

    # check para 3.4.1
    if width_eff / 2 < fd_width / 6:
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


def capacity_salgado_2008(sl, fd, h_l=0, h_b=0, vertical_load=1, verbose=0, **kwargs):
    """
    calculates the capacity according to
     THe Engineering of Foundations textbook by Salgado

     ISBN: 0072500581

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
    if not kwargs.get("disable_requires", False):
        models.check_required(sl, ["phi_r", "cohesion", "unit_dry_weight"])
        models.check_required(fd, ["length", "width", "depth"])
    v_l = kwargs.get("loc_v_l", fd.length / 2)  # given in fd coordinates
    v_b = kwargs.get("loc_v_b", fd.width / 2)
    h_eff_b = kwargs.get("h_eff_b", 0)
    h_eff_l = kwargs.get("h_eff_l", 0)
    # TODO: implement phi as fn of sigma_m

    if fd.length > fd.width:  # TODO: deal with plane strain
        fd_length = fd.length
        fd_width = fd.width
        loc_v_l = v_l
        loc_v_b = v_b
        loc_h_b = h_b
        loc_h_l = h_l
    else:
        fd_length = fd.width
        fd_width = fd.length
        loc_v_l = v_b
        loc_v_b = v_l
        loc_h_b = h_l
        loc_h_l = h_b

    ecc_b = loc_h_b * h_eff_b / vertical_load
    ecc_l = loc_h_l * h_eff_l / vertical_load

    width_eff = min(fd_width, 2 * (loc_v_b + ecc_b), 2 * (fd_width - loc_v_b - ecc_b))
    length_eff = min(fd_length, 2 * (loc_v_l + ecc_l), 2 * (fd_length - loc_v_l - ecc_l))

    # check para 3.4.1
    if width_eff / 2 < fd_width / 6:
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
    d_over_b = min(1, fd.depth / width_eff)  # limit to 1
    d_q = 1 + 2 * np.tan(sl.phi_r) * (1 - np.sin(sl.phi_r)) ** 2 * d_over_b
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


def size_footing_for_capacity(sl, vertical_load, fos=1.0, length_to_width=1.0, verbose=0, **kwargs):
    """
    Determine the size of a footing given an aspect ratio and a load

    :param sl: Soil object
    :param vertical_load: The applied load to the foundation
    :param fos: The target factor of safety
    :param length_to_width: The desired length to width ratio of the foundation
    :param verbose: verbosity
    :return: a Foundation object
    """
    method = kwargs.get("method", 'vesic')
    depth_to_width = kwargs.get("depth_to_width", 0)
    depth = kwargs.get("depth", 0)
    use_depth_to_width = 0
    if not depth:
        use_depth_to_width = 1

    # Find approximate size
    fd = models.RaftFoundation()
    fd.width = .5  # start with B=1.0m
    for i in range(50):
        fd.length = length_to_width * fd.width
        if use_depth_to_width:
            fd.depth = depth_to_width * fd.width
        else:
            fd.depth = depth
        capacity_method_selector(sl, fd, method)
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
        capacity_method_selector(sl, fd, method)
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
            capacity_method_selector(sl, fd, method)
            break
        if fs == len(fs_array) - 1:
            DesignError("No suitable foundation sizes could be determined!")

    return fd


def calc_crit_span(sl, fd, vertical_load, ip_axis='length', verbose=0, **kwargs):
    """
    Determine the size of a footing given an aspect ratio and a load

    :param sl: Soil object
    :param vertical_load: The applied load to the foundation
    :param fos: The target factor of safety
    :param length_to_width: The desired length to width ratio of the foundation
    :param verbose: verbosity
    :return: a Foundation object
    """
    method = kwargs.get("method", 'vesic')

    # Find approximate size
    new_fd = models.RaftFoundation()
    new_fd.width = fd.width
    new_fd.depth = fd.depth
    new_fd.length = fd.length
    prev_ub_len = fd.length
    q_ult = capacity_method_selector(sl, new_fd, method, verbose=max(0, verbose-1))
    init_fos = (q_ult * fd.area) / vertical_load
    if init_fos < 1.0:
        raise ValueError
    prev_lb_len = 0  # should be FOS lower than 1.
    l_ip = getattr(fd, ip_axis)
    est_len = l_ip / init_fos
    prev_q = q_ult
    for i in range(50):
        setattr(new_fd, ip_axis, est_len)
        q = capacity_method_selector(sl, new_fd, method, verbose=max(0, verbose-1))
        curr_fos = (q * new_fd.area) / vertical_load
        if np.isclose(curr_fos, 1.0, rtol=0.01):
            break
        elif curr_fos < 1:
            prev_lb_len = est_len
            est_len = (prev_ub_len + est_len) / 2
        else:
            prev_ub_len = est_len
            est_len = (prev_lb_len + est_len) / 2
    if i == 49:
        raise ValueError(init_fos, curr_fos, est_len, prev_lb_len, prev_ub_len)
    return est_len


def capacity_method_selector(sl, fd, method, **kwargs):
    """
    Calculates the bearing capacity of a foundation on soil using the specified method.
    :param sl: Soil Object
    :param fd: Foundation Object
    :param method: Method
    :param kwargs:
    :return:
    """

    if method == 'vesic':
        return capacity_vesic_1975(sl, fd, **kwargs)
    if method == 'vesics':
        return capacity_vesic_1975(sl, fd, **kwargs)
    elif method == 'nzs':
        return capacity_nzs_vm4_2011(sl, fd, **kwargs)
    elif method == 'terzaghi':
        return capacity_terzaghi_1943(sl, fd, **kwargs)
    elif method == 'hansen':
        return capacity_hansen_1970(sl, fd, **kwargs)
    elif method == 'meyerhoff':
        return capacity_meyerhof_1963(sl, fd, **kwargs)
    elif method == 'salgado':
        return capacity_salgado_2008(sl, fd, **kwargs)
    else:
        raise ValueError("method must be 'vesic', 'nzs', 'terzaghi', ")


available_methods = {
    "vesic": capacity_vesics_1975,
    "nzs": capacity_nzs_vm4_2011,
    "terzaghi": capacity_terzaghi_1943,
    "hansen": capacity_hansen_1970,
    "meyerhoff": capacity_meyerhof_1963,
    "salgado": capacity_salgado_2008
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


def capacity_meyerhof_and_hanna_1978(sl_0, sl_1, h0, fd, gwl=1e6, verbose=0):
    """
    Calculates the two-layered foundation capacity according Meyerhof and Hanna (1978)

    :param sl_0: Top Soil object
    :param sl_1: Base Soil object
    :param h0: Height of top soil layer
    :param fd: Foundation object
    :param wtl: water table level
    :param verbose: verbosity
    :return: ultimate bearing stress
    """
    sp = sm.SoilProfile()
    sp.add_layer(0, sl_0)
    sp.add_layer(h0, sl_1)
    sp.gwl = gwl
    return capacity_sp_meyerhof_and_hanna_1978(sp, fd)


def capacity_sp_meyerhof_and_hanna_1978(sp, fd, verbose=0):
    """
    Calculates the two-layered foundation capacity according Meyerhof and Hanna (1978)

    :param sp: Soil profile object
    :param fd: Foundation object
    :param wtl: water table level
    :param verbose: verbosity
    :return: ultimate bearing stress
    """
    assert isinstance(sp, sm.SoilProfile)
    if fd.length > fd.width:  # TODO: deal with plane strain
        fd_length = fd.length
        fd_width = fd.width
    else:
        fd_length = fd.width
        fd_width = fd.length

    sl_0 = sp.layer(1)
    sl_1 = sp.layer(2)
    h0 = sp.layer_depth(2)
    gwl = sp.gwl

    sl_0.nq_factor_0 = (
    (np.tan(np.pi / 4 + np.deg2rad(sl_0.phi / 2))) ** 2 * np.exp(np.pi * np.tan(np.deg2rad(sl_0.phi))))
    if sl_0.phi == 0:
        sl_0.nc_factor_0 = 5.14
    else:
        sl_0.nc_factor_0 = (sl_0.nq_factor_0 - 1) / np.tan(np.deg2rad(sl_0.phi))
    sl_0.ng_factor_0 = (sl_0.nq_factor_0 - 1) * np.tan(1.4 * np.deg2rad(sl_0.phi))

    sl_1.nq_factor_1 = (
        (np.tan(np.pi / 4 + np.deg2rad(sl_1.phi / 2))) ** 2 * np.exp(np.pi * np.tan(np.deg2rad(sl_1.phi))))
    if sl_1.phi == 0:
        sl_1.nc_factor_1 = 5.14
    else:
        sl_1.nc_factor_1 = (sl_1.nq_factor_1 - 1) / np.tan(np.deg2rad(sl_1.phi))
    sl_1.ng_factor_1 = (sl_1.nq_factor_1 - 1) * np.tan(1.4 * np.deg2rad(sl_1.phi))

    if verbose:
        log("Nc: ", sl_1.nc_factor_1)
        log("Nq: ", sl_1.nq_factor_1)
        log("Ng: ", sl_1.ng_factor_1)

    sl_0.kp_0 = (np.tan(np.pi / 4 + np.deg2rad(sl_0.phi / 2))) ** 2
    sl_1.kp_1 = (np.tan(np.pi / 4 + np.deg2rad(sl_1.phi / 2))) ** 2

    # shape factors
    if sl_0.phi >= 10:
        sl_0.s_c_0 = 1 + 0.2 * sl_0.kp_0 * (fd_width / fd_length)
        sl_0.s_q_0 = 1.0 + 0.1 * sl_0.kp_0 * (fd_width / fd_length)
    else:
        sl_0.s_c_0 = 1 + 0.2 * (fd_width / fd_length)
        sl_0.s_q_0 = 1.0
    sl_0.s_g_0 = sl_0.s_q_0

    if sl_1.phi >= 10:
        sl_1.s_c_1 = 1 + 0.2 * sl_1.kp_1 * (fd_width / fd_length)
        sl_1.s_q_1 = 1.0 + 0.1 * sl_1.kp_1 * (fd_width / fd_length)
    else:
        sl_1.s_c_1 = 1 + 0.2 * (fd_width / fd_length)
        sl_1.s_q_1 = 1.0
    sl_1.s_g_1 = sl_1.s_q_1

    # Note: this method explicitly accounts for the foundation depth, so there are no depth factors
    # TODO: inclination factors, see doi.org/10.1139/t78-060

    # Capacity
    a = 1  # assumed to  be one but can range between 1.1 and 1.27 for square footings according to Das (1999) Ch 4
    s = 1

    r = 1 + (fd_width / fd_length)

    # put the same things before that condition
    # effective weight not in the soil object

    if gwl == 0:  # case 1: GWL at surface
        q_at_interface = sl_0.unit_bouy_weight * h0
        unit_eff_weight_0_at_fd_depth = sl_0.unit_bouy_weight
        unit_eff_weight_0_at_interface = sl_0.unit_bouy_weight
        unit_eff_weight_1_below_foundation = sl_1.unit_bouy_weight

    elif 0 < gwl <= fd.depth:  # Case 2: GWL at between foundation depth and surface
        q_at_interface = (sl_0.unit_dry_weight * gwl) + (sl_0.unit_bouy_weight * (h0 - gwl))
        q_d = (sl_0.unit_dry_weight * gwl) + (sl_0.unit_bouy_weight * (fd.depth - gwl))
        unit_eff_weight_0_at_fd_depth = q_d / fd.depth
        unit_eff_weight_0_at_interface = sl_0.unit_bouy_weight
        unit_eff_weight_1_below_foundation = sl_1.unit_bouy_weight

    elif fd.depth < gwl <= fd_width + fd.depth:
        if gwl < h0:  # Case 3: GWL at between foundation depth and foundation depth plus width, and GWL < layer 1 depth

            average_unit_bouy_weight = sl_0.unit_bouy_weight + (
            ((gwl - fd.depth) / fd_width) * (sl_0.unit_dry_weight - sl_0.unit_bouy_weight))

            q_at_interface = (sl_0.unit_dry_weight * gwl) + (sl_0.unit_bouy_weight * (h0 - gwl))
            unit_eff_weight_0_at_fd_depth = sl_0.unit_dry_weight
            unit_eff_weight_0_at_interface = average_unit_bouy_weight
            unit_eff_weight_1_below_foundation = sl_1.unit_bouy_weight

        else:  # Case 4: GWL at between foundation depth and foundation depth plus width, and GWL > layer 1 depth
            average_unit_bouy_weight = sl_1.unit_bouy_weight + (
            ((gwl - h0) / fd_width) * (sl_1.unit_dry_weight - sl_1.unit_bouy_weight))

            q_at_interface = sl_0.unit_dry_weight * h0
            unit_eff_weight_0_at_fd_depth = sl_0.unit_dry_weight
            unit_eff_weight_0_at_interface = sl_0.unit_dry_weight
            unit_eff_weight_1_below_foundation = average_unit_bouy_weight

    elif gwl > fd.depth + fd_width:  # Case 5: GWL beyond foundation depth plus width
        q_at_interface = sl_0.unit_dry_weight * h0
        unit_eff_weight_0_at_fd_depth = sl_0.unit_dry_weight
        unit_eff_weight_0_at_interface = sl_0.unit_dry_weight
        unit_eff_weight_1_below_foundation = sl_1.unit_dry_weight
    else:
        raise ValueError("Could not interpret inputs")  # never reached

    # maximum value (qu <= qt)
    q_ult6 = q_at_interface - unit_eff_weight_0_at_fd_depth * fd.depth
    q_0 = (sl_0.cohesion * sl_0.nc_factor_0) + (0.5 * unit_eff_weight_0_at_interface * fd_width * sl_0.ng_factor_0)
    q_b2 = (q_at_interface * sl_1.nq_factor_1 * sl_1.s_q_1)
    q_1 = (sl_1.cohesion * sl_1.nc_factor_1) + (0.5 * unit_eff_weight_1_below_foundation * fd_width * sl_1.ng_factor_1)
    q_b3 = (unit_eff_weight_1_below_foundation * fd_width * sl_1.ng_factor_1 * sl_1.s_g_1 / 2)
    q_ult5 = r * (unit_eff_weight_0_at_interface * ((h0 - fd.depth) ** 2)) * (1 + (2 * fd.depth / (h0 - fd.depth))) * (
        np.tan(np.deg2rad(sl_0.phi)) / fd_width) * s
    q_t2 = (unit_eff_weight_0_at_fd_depth * fd.depth * sl_0.nq_factor_0 * sl_0.s_q_0)
    q_t3 = (unit_eff_weight_0_at_interface * fd_width * sl_0.ng_factor_0 * sl_0.s_g_0 / 2)

    # qb
    q_b1 = (sl_1.cohesion * sl_1.nc_factor_1 * sl_1.s_c_1)
    q_b = q_b1 + q_b2 + q_b3

    q1_q0 = q_1 / q_0

    # calculate the ca factor
    # if sl_0.cohesion == 0:
    #     c1_c0 = 0
    # else:
    #     c1_c0 = sl_1.cohesion / sl_0.cohesion
    x = np.array([0.000, 0.082, 0.206, 0.298, 0.404, 0.509, 0.598, 0.685, 0.772])
    y = np.array([0.627, 0.700, 0.794, 0.855, 0.912, 0.948, 0.968, 0.983, 0.997])

    # raise Warning("ca should be interpolated using q1/q2 not cohesion, see Figure 4 in MH1978")
    ca_c0 = np.interp(q1_q0, x, y)
    ca = ca_c0 * sl_0.cohesion

    # ks
    x_0 = np.array([0, 20.08, 22.42, 25.08, 27.58, 30.08, 32.58, 34.92, 37.83, 40.00, 42.67, 45.00, 47.00, 49.75])
    y_0 = np.array([0.93, 0.93, 0.93, 0.93, 1.01, 1.17, 1.32, 1.56, 1.87, 2.26, 2.72, 3.35, 3.81, 4.82])
    x_2 = np.array([0, 20.08, 22.50, 25.08, 27.58, 30.08, 32.50, 35.00, 37.67, 40.17, 42.67, 45.00, 47.50, 50.00])
    y_2 = np.array([1.55, 1.55, 1.71, 1.86, 2.10, 2.33, 2.72, 3.11, 3.81, 4.43, 5.28, 6.14, 7.46, 9.24])
    x_4 = np.array([0, 20.00, 22.51, 25.10, 27.69, 30.11, 32.45, 35.04, 37.88, 40.14, 42.65, 45.07, 47.33, 50.08])
    y_4 = np.array([2.49, 2.49, 2.64, 2.87, 3.34, 3.81, 4.43, 5.20, 6.29, 7.38, 9.01, 11.11, 14.29, 19.34])
    x_10 = np.array([0, 20.00, 22.50, 25.08, 28.00, 30.00, 32.50, 34.92, 37.50, 40.17, 42.42, 45.00, 47.17, 50.08])
    y_10 = np.array([3.27, 3.27, 3.74, 4.44, 5.37, 6.07, 7.16, 8.33, 10.04, 12.30, 15.95, 21.17, 27.47, 40.00])
    x_int = sl_0.phi

    if sl_0.phi < 1:
        fd.ks = 0
    else:

        if q1_q0 == 0:
            fd.ks = np.interp(x_int, x_0, y_0)

        elif q1_q0 == 0.2:
            fd.ks = np.interp(x_int, x_2, y_2)

        elif q1_q0 == 0.4:
            fd.ks = np.interp(x_int, x_4, y_4)

        elif q1_q0 == 1.0:
            fd.ks = np.interp(x_int, x_10, y_10)

        elif 0 < q1_q0 < 0.2:
            ks_1 = np.interp(x_int, x_0, y_0)
            ks_2 = np.interp(x_int, x_2, y_2)
            fd.ks = (((ks_2 - ks_1) * q1_q0) / 0.2) + ks_1

        elif 0.2 < q1_q0 < 0.4:
            ks_1 = np.interp(x_int, x_2, y_2)
            ks_2 = np.interp(x_int, x_4, y_4)
            fd.ks = (((ks_2 - ks_1) * (q1_q0 - 0.2)) / 0.2) + ks_1

        elif 0.4 < q1_q0 < 1.0:
            ks_1 = np.interp(x_int, x_4, y_4)
            ks_2 = np.interp(x_int, x_10, y_10)
            fd.ks = (((ks_2 - ks_1) * (q1_q0 - 0.4)) / 0.6) + ks_1
        else:
            raise DesignError(
                "Cannot compute 'ks', bearing ratio out-of-range (q1_q0 = %.3f) required: 0-1." % q1_q0)

    # qu
    q_ult4 = (r * (2 * ca * (h0 - fd.depth) / fd_width) * a)
    q_ult5_ks = q_ult5 * fd.ks
    q_ult = q_b + q_ult4 + q_ult5_ks - q_ult6

    q_t1 = (sl_0.cohesion * sl_0.nc_factor_0 * sl_0.s_c_0)
    q_t = q_t1 + q_t2 + q_t3

    if q_ult > q_t:
        if h0 > fd_width/2:
            fd.q_ult = q_t

        else:
            vert_eff_stress_interface = sp.get_v_eff_stress_at_depth(h0)
            vert_eff_stress_lowest = sp.get_v_eff_stress_at_depth(fd_width+fd.depth)
            average_eff_stress = (vert_eff_stress_interface + vert_eff_stress_lowest) / 2

            c_2_eff = sl_1.cohesion + average_eff_stress * np.tan(np.radians(sl_1.phi))

            if sl_0.cohesion > c_2_eff:
                fd.q_ult = q_t

            else:
                # vd = {}
                # vd[1] =[1, 1, 1, 1, 1]
                # vd[0.667] = [1, 1.033, 1.064, 1.088, 1.109]
                # vd[0.5] = [1, 1.056, 1.107, 1.152, 1.193]
                # vd[0.333] = [1, 1.088, 1.167, 1.241, 1.311]
                # vd[0.25] = [1, 1.107, 1.208, 1.302, 1.389]
                # vd[0.2] = [1, 1.121, 1.235, 1.342, 1.444]
                # vd[0.1] = [1, 1.154, 1.302, 1.446, 1.584]

                h_over_b = (h0 - fd.depth) / fd_width
                c1_over_c2 = sl_0.cohesion/c_2_eff

                c_1_over_c_2 = [0.1, 0.2, 0.25, 0.333, 0.5, 0.667, 1.]
                m_1 = [1.584, 1.444, 1.389, 1.311, 1.193, 1.109, 1.]
                m_125 = [1.446, 1.342, 1.302, 1.241, 1.152, 1.088, 1.]
                m_167 = [1.302, 1.235, 1.208, 1.167, 1.107, 1.064, 1.]
                m_25 = [1.154, 1.121, 1.107, 1.088, 1.056, 1.033, 1.]
                m_5 = [1, 1, 1, 1, 1, 1, 1]

                if h_over_b == 0.1:
                    m = np.interp(c1_over_c2, c_1_over_c_2, m_1)
                elif h_over_b == 0.125:
                    m = np.interp(c1_over_c2, c_1_over_c_2, m_125)
                elif h_over_b == 0.167:
                    m = np.interp(c1_over_c2, c_1_over_c_2, m_167)
                elif h_over_b == 0.250:
                    m = np.interp(c1_over_c2, c_1_over_c_2, m_25)
                elif h_over_b >= 0.5:
                    m = np.interp(c1_over_c2, c_1_over_c_2, m_5)
                elif 0.1 < h_over_b < 0.125:
                    m_a = np.interp(c1_over_c2, c_1_over_c_2, m_1)
                    m_b = np.interp(c1_over_c2, c_1_over_c_2, m_125)
                    m = np.interp(h_over_b, [0.1,0.125], [m_a,m_b])
                elif 0.125 < h_over_b < 0.167:
                    m_a = np.interp(c1_over_c2, c_1_over_c_2, m_125)
                    m_b = np.interp(c1_over_c2, c_1_over_c_2, m_167)
                    m = np.interp(h_over_b, [0.125, 0.167], [m_a, m_b])
                elif 0.167 < h_over_b < 0.25:
                    m_a = np.interp(c1_over_c2, c_1_over_c_2, m_167)
                    m_b = np.interp(c1_over_c2, c_1_over_c_2, m_25)
                    m = np.interp(h_over_b, [0.167, 0.250], [m_a, m_b])
                elif 0.25 < h_over_b < 0.5:
                    m_a = np.interp(c1_over_c2, c_1_over_c_2, m_25)
                    m_b = np.interp(c1_over_c2, c_1_over_c_2, m_5)
                    m = np.interp(h_over_b, [0.250, 0.500], [m_a, m_b])

                fd.q_ult = (sl_0.cohesion * m * sl_0.nc_factor_0) + (unit_eff_weight_0_at_fd_depth * fd.depth)

    else:
        fd.q_ult = q_ult

    return fd.q_ult




def calc_norm_moment_bs_via_butterfield_et_al_1994(n, v, n_max, b, qv_max, qm_max):
    c = 0.22
    t = qv_max
    s = qm_max * b
    m = 0.5 * (2 * c * s * v / t - s ** 2 * np.sqrt((4 * c ** 2 * v ** 2 / (s ** 2 * t ** 2) -
        4 * (-n ** 4 / n_max ** 2 + 2 * n ** 3 / n_max - n ** 2 + v ** 2 / t ** 2) / s ** 2)))

    m = c * s * v / t - s * (np.sqrt((c ** 2 - 1) * v ** 2 / t ** 2 + n ** 2 * (n - n_max) ** 2 / n_max ** 2)) / t
    # m = s * np.sqrt(c ** 2 * v ** 2 + n ** 4 * t ** 2 - 2 * n ** 3 * t ** 2 - v ** 2) / t + c * s * v / t
    return m


def calc_norm_shear_bs_via_butterfield_et_al_1994(n, r, qv_max, qm_max):
    c = 0.22
    t = qv_max
    s = qm_max
    v2 = (n - 1) * n * s * t / (np.sqrt(2 * c * r * s * t - r ** 2 * t ** 2 - s ** 2))
    return v2

def calc_norm_shear_bs(n, m):
    c = 0.22
    u = (c ** 2 - 1) * m ** 2 + (n - 1) ** 2 * n ** 2


def run():
    import matplotlib.pyplot as plt
    n = 0.5
    n_max = 1
    b = 1
    qv_max = 0.52
    qm_max = 0.35
    v = np.linspace(0.1, 0.14, 20)
    # m = calc_norm_moment_bs_via_butterfield_et_al_1994(n, v, n_max, b, qv_max, qm_max)
    # for i in range(len(m)):
    #     print(i, v[i], m[i])
    # plt.plot(v / n_max, m / (n_max))
    v = 0
    n = np.linspace(0.1, 0.9, 10)
    n = 0.5
    r = np.logspace(-1, 2, 5)
    qv = calc_norm_shear_bs_via_butterfield_et_al_1994(n, r, qv_max, qm_max)
    print(qv)
    qm = r * qv
    plt.plot(qv, qm)
    # m = calc_norm_moment_bs_via_butterfield_et_al_1994(n, v, n_max, b, qv_max, qm_max)
    # plt.plot(n / n_max, m / (b * n_max))
    plt.show()




if __name__ == '__main__':
    run()