import numpy as np

from geofound.output import log
from geofound.exceptions import DesignError
from geofound import models
import sfsimodels as sm


def calc_shape_factor_coh_salgado_et_al_2004(fd):
    """
    Salgado (2008) Eq. 10.22 and Table 10-3

    Parameters
    ----------
    fd: sm.Foundation

    Returns
    -------

    """
    if fd.length > fd.width:
        fd_length = fd.length
        fd_width = fd.width
    else:
        fd_length = fd.width
        fd_width = fd.length
    b_o_l = fd_width / fd_length
    bols = [0.2, 0.25, 0.33, 0.5, 1]
    c1s = [0.190, 0.172, 0.159, 0.156, 0.125]
    c2s = [0.090, 0.11, 0.137, 0.173, 0.219]
    c1 = np.interp(b_o_l, bols, c1s)
    c2 = np.interp(b_o_l, bols, c2s)
    return 1 + c1 * b_o_l + c2 * np.sqrt(fd.depth / fd_width)


def calc_shape_factor_coh_gourvenec_2007(fd):
    if fd.length > fd.width:
        fd_length = fd.length
        fd_width = fd.width
    else:
        fd_length = fd.width
        fd_width = fd.length
    b_o_l = fd_width / fd_length
    return 1 + 0.214 * b_o_l - 0.067 * b_o_l ** 2


def calc_phi_peak_bolton_1986(phi_c_txc, p_eff, d_r, k=2, q=10, r=1, p_atm=101.0e3):
    """

    :param phi_c_txc: float
        Constant volume friction angle in triaxial compression [degrees]
    :param p_eff:
    :param sl:
    :param k: k=1 for plane strain, =2 for triaxial compression
    :param q:
    :param r:
    :param p_atm:
    :return:
    """
    i_r = d_r * (q - np.log(100 * p_eff / p_atm)) - r  # From Salgado 2000 (adapted to give consistent units)
    # i_r = sl.relative_density * (q - np.log(p_eff)) - r
    a_psi = (5 - 2 * (k - 1))
    return phi_c_txc + a_psi * i_r


def calc_phi_peak_fd_salgado_2008(phi_c_txc, p_eff, d_r, l_o_b=None, fd=None, q=10, r=1, p_atm=101.0e3,
                                  phi_c_ps_diff=0.0):
    """

    :param phi_c_txc: float
        Constant volume friction angle in triaxial compression [degrees]
    :param p_eff:
    :param sl:
    :param k: k=1 for plane strain, =2 for triaxial compression
    :param q:
    :param r:
    :param p_atm:
    :return:
    """
    if l_o_b is None:
        l_o_b = max(fd.length / fd.width, fd.width / fd.length)
    i_r = d_r * (q - np.log(100 * p_eff / p_atm)) - r  # From Salgado 2000 (adapted to give consistent units)
    # i_r = sl.relative_density * (q - np.log(p_eff)) - r
    a_psi = np.clip(1. / 3 * (l_o_b + 8), 3, 5)  # This is originally from Perkins and Madson (2000)
    phi_c = phi_c_txc + (a_psi - 3) / 2 * phi_c_ps_diff
    return phi_c + a_psi * i_r


def calc_phi_bc_strip_loukidis_2019(phi_c_txc, sl, fd, ip_axis_2d, p_atm=101.0e3, gwl=1e6):
    """
    Frictional angle for calculating bearing capacity of a strip footing

    Parameters
    ----------
    phi_c_txc: float
        Critical state friction angle under triaxial compression
    sl
    fd
    ip_axis_2d
    p_atm
    gwl

    Returns
    -------

    """
    if ip_axis_2d == 'length':
        fd_width = fd.length
    elif ip_axis_2d == 'width':
        fd_width = fd.width
    else:
        raise ValueError(f'ip_axis_2d must be either: "width", or "length" not {ip_axis_2d}')
    if gwl == 0:
        unit_weight = sl.unit_bouy_weight
    elif gwl >= fd.depth + fd_width:
        unit_weight = sl.unit_dry_weight
    elif 0 < gwl < fd.depth:
        unit_weight = sl.unit_bouy_weight
    elif fd.depth <= gwl <= fd.depth + fd_width:
        average_unit_bouy_weight = sl.unit_bouy_weight + (
                ((gwl - fd.depth) / fd_width) * (sl.unit_dry_weight - sl.unit_bouy_weight))
        unit_weight = average_unit_bouy_weight
    return phi_c_txc + ((17.6 * sl.relative_density - 8.8) - 2.44 * np.log(fd_width * unit_weight / p_atm))


def calc_m_eff_bc_via_loukidis_and_salgado_2006(sl, fd, ip_axis_2d=None, p_atm=101.0e3, gwl=1e6, ob=0.0):

    if ip_axis_2d is None:
        if fd.length > fd.width:
            fd_length = fd.length
            fd_width = fd.width
        else:
            fd_length = fd.width
            fd_width = fd.length
        b_o_l = fd_width / fd_length
    elif ip_axis_2d == 'length':
        fd_width = fd.length
        b_o_l = 1
    elif ip_axis_2d == 'width':
        fd_width = fd.width
        b_o_l = 1
    else:
        raise ValueError(f'ip_axis_2d must be either: None, "width", or "length" not {ip_axis_2d}')
    assert 0.01 < sl.unit_weight / p_atm < 10, (sl.unit_weight, p_atm)
    if gwl == 0:
        unit_weight = sl.unit_bouy_weight
    elif gwl >= fd.depth + fd_width:
        unit_weight = sl.unit_dry_weight
    elif 0 < gwl < fd.depth:
        unit_weight = sl.unit_bouy_weight
    elif fd.depth <= gwl <= fd.depth + fd_width:
        average_unit_bouy_weight = sl.unit_bouy_weight + (
                ((gwl - fd.depth) / fd_width) * (sl.unit_dry_weight - sl.unit_bouy_weight))
        unit_weight = average_unit_bouy_weight
    sigma_meff = 20 * p_atm * ((unit_weight * fd_width + ob) / p_atm) ** 0.7 * (1 - 0.32 * b_o_l)
    return sigma_meff


def calc_m_eff_bc_via_debeer_1965(sl, fd, q_ult, ip_axis_2d=None, gwl=1e6):
    """
    Calculation of average confining stress to compute peak friction angle for ultimate bearing capacity

    # Equation taken from Salgado 2008 page 442

    Parameters
    ----------
    sl
    fd
    q_ult
    ip_axis_2d
    gwl

    Returns
    -------

    """
    if ip_axis_2d is None:
        if fd.length > fd.width:
            fd_width = fd.width
        else:
            fd_width = fd.length
    elif ip_axis_2d == 'length':
        fd_width = fd.length
    elif ip_axis_2d == 'width':
        fd_width = fd.width
    else:
        raise ValueError(f'ip_axis_2d must be either: None, "width", or "length" not {ip_axis_2d}')
    if gwl == 0:
        q_d = sl.unit_eff_weight * fd.depth
    elif 0 < gwl < fd.depth:
        q_d = (sl.unit_dry_weight * gwl) + (sl.unit_bouy_weight * (fd.depth - gwl))

    elif fd.depth <= gwl <= fd.depth + fd_width:
        sl.average_unit_bouy_weight = sl.unit_bouy_weight + (
                ((gwl - fd.depth) / fd_width) * (sl.unit_dry_weight - sl.unit_bouy_weight))
        q_d = sl.unit_dry_weight * fd.depth
    sigma_meff = 1. / ((1 + np.tan(sl.phi_r) ** 2) * (1 + np.sin(sl.phi_r))) * (q_ult + 3 * q_d) / 4
    return sigma_meff


def calc_m_eff_bc_via_perkins_and_madson_2000(fd, q_demand, ip_axis_2d=None, min_lim=True):
    # for bearing capacity
    if ip_axis_2d is None:
        if fd.length > fd.width:
            fd_length = fd.length
            fd_width = fd.width
        else:
            fd_length = fd.width
            fd_width = fd.length
        l_o_b = fd_length / fd_width
    elif ip_axis_2d == 'length':
        fd_width = fd.length
        l_o_b = 7  # pg 526 PM2000
    elif ip_axis_2d == 'width':
        fd_width = fd.width
        l_o_b = 7
    else:
        raise ValueError(f'ip_axis_2d must be either: None, "width", or "length" not {ip_axis_2d}')
    a = 1. / 6 * (0.52 - 0.04 * l_o_b)
    if min_lim:
        sigma_meff = max(1. / 6 * (0.52 - 0.04 * l_o_b) * q_demand, q_demand / 25)
    else:
        sigma_meff = 1. / 6 * (0.52 - 0.04 * l_o_b) * q_demand
    return sigma_meff


def capacity_vesic_1975(sl, fd, slope=0, base_tilt=0, ip_axis_2d=None, verbose=0,
                        gwl=1e6, **kwargs):
    """
    Calculates the foundation capacity according Vesics(1975)
    #Gunaratne, Manjriker. 2006. "Spread Footings: Analysis and Design."
    Ref: http://geo.cv.nctu.edu.tw/foundation/download/
                            BearingCapacityOfFoundations.pdf

    Note: Phi should be from plane strain (i.e. higher than if rectangular)

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
    if 'axis_inf' in kwargs:
        raise ValueError('use ip_axis_2d instead of axis_inf')

    hload_l = kwargs.setdefault('hload_length', 0)
    hload_b = kwargs.setdefault('hload_width', 0)
    nload = kwargs.setdefault('nload', 1)
    horizontal_load = np.sqrt(hload_l ** 2 + hload_b ** 2)
    # eccentricity due to offset vertical load, hload at height and direct moment.
    e_length = kwargs.setdefault('e_length', 0)
    e_width = kwargs.setdefault('e_width', 0)

    ob = kwargs.get('ob', 0.0)

    area_foundation = fd.length * fd.width
    temp_fd_length = fd.length
    temp_fd_width = fd.width
    if ip_axis_2d is not None:
        if ip_axis_2d == 'width':
            area_foundation = temp_fd_width
            temp_fd_length = temp_fd_width * 100
        elif ip_axis_2d == 'length':
            area_foundation = temp_fd_length
            temp_fd_width = temp_fd_length * 100
        else:
            raise ValueError(f"ip_axis_2d must be either 'width' or 'length', not {ip_axis_2d}")

    temp_width_eff = temp_fd_width - 2 * abs(e_width)
    temp_length_eff = temp_fd_length - 2 * abs(e_length)
    if temp_fd_length > temp_fd_width:
        fd_length = temp_fd_length
        fd_width = temp_fd_width
        rev = 0
    else:
        fd_length = temp_fd_width
        fd_width = temp_fd_length
        rev = 1

    width_for_stress = min(temp_width_eff, temp_length_eff)
    c_a = 0.6 - 1.0 * sl.cohesion

    fd.nq_factor = ((np.tan(np.pi / 4 + sl.phi_r / 2)) ** 2 * np.exp(np.pi * np.tan(sl.phi_r)))
    if sl.phi_r == 0:
        fd.nc_factor = 5.14
    else:
        fd.nc_factor = (fd.nq_factor - 1) / np.tan(sl.phi_r)
    fd.ng_factor = 2.0 * (fd.nq_factor + 1) * np.tan(sl.phi_r)

    # Note: according to Salgado (2008) page 437 Vesic (1973) uses true dimensions for shape and depth factors
    # shape factors:
    if fd_length / fd_width > 10:
        s_c = 1.0
        s_q = 1.0
        s_g = 1.0
    else:
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
    if fd_length / fd_width > 10:
        m = 2.0
    else:
        m_b = (2.0 + fd_width / fd_length) / (1 + fd_width / fd_length)
        m_l = (2.0 + fd_length / fd_width) / (1 + fd_length / fd_width)
        m = np.sqrt(m_b ** 2 + m_l ** 2)

    if sl.phi_r == 0:
        i_q = 1.0
        i_g = 1.0
    else:
        i_q = (1.0 - horizontal_load / (nload + area_foundation *
                                        c_a / np.tan(sl.phi_r))) ** m
        i_g = (1.0 - horizontal_load / (nload + area_foundation *
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
    elif 0 < gwl < fd.depth:
        q_d = (sl.unit_dry_weight * gwl) + (sl.unit_bouy_weight * (fd.depth - gwl))
        unit_weight = sl.unit_bouy_weight
    elif fd.depth <= gwl <= fd.depth + fd_width:
        sl.average_unit_bouy_weight = sl.unit_bouy_weight + (
                ((gwl - fd.depth) / fd_width) * (sl.unit_dry_weight - sl.unit_bouy_weight))
        q_d = sl.unit_dry_weight * fd.depth
        unit_weight = sl.average_unit_bouy_weight
    elif gwl > fd.depth + fd_width:
        q_d = sl.unit_dry_weight * fd.depth
        unit_weight = sl.unit_dry_weight
    else:
        raise ValueError(f'gwl must be zero or positive, not {gwl}')
    q_d += ob
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
                0.5 * width_for_stress * unit_weight *
                fd.ng_factor * s_g * d_g * i_g * g_g * b_g)

    fd.width_eff = temp_width_eff
    fd.length_eff = temp_length_eff
    if ip_axis_2d is not None:
        setattr(fd, f'{fd.ip_axis}_eff', 1)
    fd.area_eff = fd.width_eff * fd.length_eff
    fd.n_ult = fd.q_ult * fd.area_eff

    if verbose:
        log("qult: ", fd.q_ult)
    return fd.q_ult


def capacity_terzaghi_1943(sl, fd, round_footing=False, ip_axis_2d=None, verbose=0, **kwargs):
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
    if 'axis_inf' in kwargs:
        raise ValueError('use ip_axis_2d instead of axis_inf')
    if ip_axis_2d is None:
        if fd.length > fd.width:
            fd_length = fd.length
            fd_width = fd.width
        else:
            fd_length = fd.width
            fd_width = fd.length
    elif ip_axis_2d == 'length':
        fd_width = fd.length
        fd_length = None
    elif ip_axis_2d == 'width':
        fd_width = fd.width
        fd_length = None
    else:
        raise ValueError(f'ip_axis_2d must be either: None, "width", or "length" not {ip_axis_2d}')

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


def capacity_brinch_hansen_1970(sl, fd, gwl=1e6, slope=0, base_tilt=0, ip_axis_2d=None, verbose=0,
                                **kwargs):
    """
    Calculates the foundation capacity according Hansen (1970)
    Ref: http://bestengineeringprojects.com/civil-projects/
    hansens-bearing-capacity-theory/

    Note: Phi should be from plane strain (i.e. higher than if rectangular)

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

    hload_l = kwargs.setdefault('hload_length', 0)
    hload_b = kwargs.setdefault('hload_width', 0)
    nload = kwargs.setdefault('nload', 1)
    horizontal_load = np.sqrt(hload_l ** 2 + hload_b ** 2)
    # eccentricity due to offset vertical load, hload at height and direct moment.
    e_length = kwargs.setdefault('e_length', 0)
    e_width = kwargs.setdefault('e_width', 0)
    ob = kwargs.get('ob', 0.0)  # overburden pressure

    if 'axis_inf' in kwargs:
        raise ValueError('use ip_axis_2d instead of axis_inf')

    area_foundation = fd.length * fd.width
    temp_fd_length = fd.length
    temp_fd_width = fd.width
    if ip_axis_2d is not None:
        if ip_axis_2d == 'width':
            area_foundation = temp_fd_width
            temp_fd_length = temp_fd_width * 100
        elif ip_axis_2d == 'length':
            area_foundation = temp_fd_length
            temp_fd_width = temp_fd_length * 100
        else:
            raise ValueError(f"ip_axis_2d must be either 'width' or 'length', not {ip_axis_2d}")

    temp_width_eff = temp_fd_width - 2 * abs(e_width)
    temp_length_eff = temp_fd_length - 2 * abs(e_length)
    if temp_length_eff > temp_width_eff:
        length_eff = temp_length_eff
        width_eff = temp_width_eff
        rev = 0
    else:
        length_eff = temp_width_eff
        width_eff = temp_length_eff
        rev = 1

    c_a = 0.6 - 1.0 * sl.cohesion

    # Note: exact solution for associated flow rule
    fd.nq_factor = ((np.tan(np.pi / 4 + sl.phi_r / 2)) ** 2 * np.exp(np.pi * np.tan(sl.phi_r)))
    if sl.phi_r == 0:
        fd.nc_factor = 5.14
    else:
        fd.nc_factor = (fd.nq_factor - 1) / np.tan(sl.phi_r)
    fd.ng_factor = 1.5 * (fd.nq_factor - 1) * np.tan(sl.phi_r)

    # Note: according to Salgado (2008) page 437 BH1970 uses effective dimensions for shape and depth factors
    # shape factors
    if length_eff / width_eff > 10:
        s_c = 1
        s_q = 1
        s_g = 1
    else:
        if sl.phi_r == 0:
            s_c = 0.2 * width_eff / length_eff
        else:
            s_c = 1.0 + fd.nq_factor / fd.nc_factor * width_eff / length_eff

        s_q = 1.0 + width_eff / length_eff * np.sin(sl.phi_r)
        s_g = 1.0 - 0.4 * width_eff / length_eff

    # depth factors:
    if fd.depth / width_eff > 1:
        k = np.arctan(fd.depth / width_eff)
    else:
        k = fd.depth / width_eff
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
                (nload + area_foundation * c_a / np.tan(sl.phi_r))) ** 5)
        i_c = i_q - (1 - i_q) / (fd.nq_factor - 1)
        i_g = ((1 - (0.7 * horizontal_load) /
                (nload + area_foundation * c_a / np.tan(sl.phi_r))) ** 5)

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
    if gwl == 0:
        q_d = sl.unit_bouy_weight * fd.depth
        unit_weight = sl.unit_bouy_weight
    elif 0 < gwl < fd.depth:
        q_d = (sl.unit_dry_weight * gwl) + (sl.unit_bouy_weight * (fd.depth - gwl))
        unit_weight = sl.unit_bouy_weight
    elif fd.depth <= gwl <= fd.depth + width_eff:
        sl.average_unit_bouy_weight = sl.unit_bouy_weight + (
                ((gwl - fd.depth) / width_eff) * (sl.unit_dry_weight - sl.unit_bouy_weight))
        q_d = sl.unit_dry_weight * fd.depth
        unit_weight = sl.average_unit_bouy_weight
    elif gwl > fd.depth + width_eff:
        q_d = sl.unit_dry_weight * fd.depth
        unit_weight = sl.unit_dry_weight
    else:
        raise DesignError(f'gwl must be positive or zero, not: {gwl}')
    q_d += ob  # add overburden

    # Capacity
    if sl.phi_r == 0:
        fd.q_ult = (sl.cohesion * fd.nc_factor *
                    (1 + s_c + d_c - i_c - g_c - b_c) + q_d)
    else:
        fd.q_ult = (sl.cohesion * fd.nc_factor *
                    s_c * d_c * i_c * g_c * b_c +
                    q_d * fd.nq_factor * s_q * d_q * i_q * g_q * b_q +
                    0.5 * width_eff * unit_weight *
                    fd.ng_factor * s_g * d_g * i_g * g_g * b_g)
    if rev:
        fd.width_eff = length_eff
        fd.length_eff = width_eff
    else:
        fd.width_eff = width_eff
        fd.length_eff = length_eff
    if ip_axis_2d is not None:
        setattr(fd, f'{fd.ip_axis}_eff', 1)
    fd.area_eff = fd.width_eff * fd.length_eff
    fd.n_ult = fd.q_ult * fd.area_eff
    if verbose:
        log("qult: ", fd.q_ult)
    return fd.q_ult


def capacity_meyerhof_1963(sl, fd, gwl=1e6, ip_axis_2d=None, verbose=0, **kwargs):
    """
    Calculates the foundation capacity according Meyerhoff (1963)
    http://www.engs-comp.com/meyerhof/index.shtml

    Phi should be appropriate for the shape (not plane strain phi, i.e. should be lower if square than strip)

    :param sl: Soil object
    :param fd: Foundation object
    :param hload_l: Horizontal load parallel to length
    :param hload_b: Horizontal load parallel to width
    :param nload: Vertical load
    :param verbose: verbosity
    :return: ultimate bearing stress
    """

    if not kwargs.get("disable_requires", False):
        models.check_required(sl, ["phi_r", "cohesion", "unit_dry_weight"])
        models.check_required(fd, ["length", "width", "depth"])
    hload_b = kwargs.setdefault('hload_length', 0)
    hload_l = kwargs.setdefault('hload_width', 0)
    nload = kwargs.setdefault('nload', 1)
    horizontal_load = np.sqrt(hload_l ** 2 + hload_b ** 2)
    # eccentricity due to offset vertical load, hload at height and direct moment.
    e_length = kwargs.setdefault('e_length', 0)
    e_width = kwargs.setdefault('e_width', 0)
    ob = kwargs.get('ob', 0.0)
    if 'axis_inf' in kwargs:
        raise ValueError('use ip_axis_2d instead of axis_inf')
    temp_fd_length = fd.length
    temp_fd_width = fd.width
    if ip_axis_2d is not None:
        if ip_axis_2d == 'width':
            temp_fd_length = temp_fd_width * 100
        elif ip_axis_2d == 'length':
            temp_fd_width = temp_fd_length * 100
        else:
            raise ValueError(f"ip_axis_2d must be either 'width' or 'length', not {ip_axis_2d}")

    temp_width_eff = temp_fd_width - 2 * abs(e_width)
    temp_length_eff = temp_fd_length - 2 * abs(e_length)
    if temp_length_eff > temp_width_eff:
        length_eff = temp_length_eff
        width_eff = temp_width_eff
        rev = 0
    else:
        length_eff = temp_width_eff
        width_eff = temp_length_eff
        rev = 1

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
    if length_eff / width_eff > 10:
        s_c = 1
        s_q = 1
    else:
        s_c = 1 + 0.2 * kp * width_eff / length_eff
        if sl.phi > 10:
            s_q = 1.0 + 0.1 * kp * width_eff / length_eff
        else:
            s_q = 1.0
    s_g = s_q

    # depth factors
    d_c = 1 + 0.2 * np.sqrt(kp) * fd.depth / width_eff
    if sl.phi > 10:
        d_q = 1 + 0.1 * np.sqrt(kp) * fd.depth / width_eff
    else:
        d_q = 1.0
    d_g = d_q

    # inclination factors:
    theta_load = np.arctan(horizontal_load / nload)
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
    elif 0 < gwl < fd.depth:
        q_d = (sl.unit_dry_weight * gwl) + (sl.unit_bouy_weight * (fd.depth - gwl))
        unit_weight = sl.unit_bouy_weight
    elif fd.depth <= gwl <= fd.depth + width_eff:
        sl.average_unit_bouy_weight = sl.unit_bouy_weight + (
                ((gwl - fd.depth) / width_eff) * (sl.unit_dry_weight - sl.unit_bouy_weight))
        q_d = sl.unit_dry_weight * fd.depth
        unit_weight = sl.average_unit_bouy_weight
    elif gwl > fd.depth + width_eff:
        q_d = sl.unit_dry_weight * fd.depth
        unit_weight = sl.unit_dry_weight
    else:
        raise DesignError(f'gwl must be positive or zero, not: {gwl}')
    q_d += ob  # add overburden

    if verbose:
        log("Nc: ", fd.nc_factor)
        log("Nq: ", fd.nq_factor)
        log("Ng: ", fd.ng_factor)
        log("s_c: %.3f" % s_c)
        log("s_q: %.3f" % s_q)
        log("s_g: %.3f" % s_g)
        log("d_c: %.3f" % d_c)
        log("d_q: %.3f" % d_q)
        log("d_g: %.3f" % d_g)
        log("i_c: %.3f" % i_c)
        log("i_q: %.3f" % i_q)
        log("i_g: %.3f" % i_g)
        log("q_d: %.3f" % q_d)

    # Capacity
    fd.q_ult = (sl.cohesion * fd.nc_factor * s_c * d_c * i_c +
                q_d * fd.nq_factor * s_q * d_q * i_q +
                0.5 * width_eff * unit_weight *
                fd.ng_factor * s_g * d_g * i_g)
    if rev:
        fd.width_eff = length_eff
        fd.length_eff = width_eff
    else:
        fd.width_eff = width_eff
        fd.length_eff = length_eff
    if ip_axis_2d is not None:
        setattr(fd, f'{fd.ip_axis}_eff', 1)
    fd.area_eff = fd.width_eff * fd.length_eff
    fd.n_ult = fd.q_ult * fd.area_eff
    if verbose:
        log("q_ult: %.3f" % fd.q_ult)
    return fd.q_ult


def capacity_nzs_vm4_2011(sl, fd, gwl=1e6, slope=0, verbose=0, save_factors=0, **kwargs):
    """
    calculates the capacity according to
     Appendix B verification method 4 of the NZ building code

    :param sl: Soil object
    :param fd: Foundation object
    :param h_loads:
        Horizontal loads dictionary 'length_dir' and 'width_dir', each containing:
            'load_ratio': the ratio of horizontal load to vertical load
                (note positive should be consistent with a eccens and right-hand-rule)
            'height': The height of the applied load
    :param e_length:
        Vertical load eccentricity in length dir
    :param e_width:
        Vertical load eccentricity in width dir
    :param slope: ground slope
    :param verbose: verbosity
    :return: ultimate bearing stress
    """
    # Need to make adjustments if sand  has DR<40% or
    # clay has liquidity indices greater than 0.7

    if not kwargs.get("disable_requires", False):
        models.check_required(sl, ["phi_r", "cohesion", "unit_dry_weight"])
        models.check_required(fd, ["length", "width", "depth"])

    hload_l = kwargs.setdefault('hload_length', 0)
    hload_b = kwargs.setdefault('hload_width', 0)
    nload = kwargs.setdefault('nload', 1)
    horizontal_load = np.sqrt(hload_l ** 2 + hload_b ** 2)
    dep_inputs = ['h_b', 'h_l', 'vertical_load', 'h_eff_b', 'h_eff_l']
    for item in dep_inputs:
        if item in kwargs:
            raise ValueError(f'Input {item} has deprecated')

    # eccentricity due to offset vertical load, hload at height and direct moment.
    e_length = kwargs.setdefault('e_length', 0)
    e_width = kwargs.setdefault('e_width', 0)

    ip_axis_2d = kwargs.get('ip_axis_2d', None)
    ob = kwargs.get('ob', 0.0)  # overburden pressure

    temp_fd_length = fd.length
    temp_fd_width = fd.width
    if ip_axis_2d is not None:
        if ip_axis_2d == 'width':
            temp_fd_length = temp_fd_width * 100
        elif ip_axis_2d == 'length':
            temp_fd_width = temp_fd_length * 100
        else:
            raise ValueError(f"ip_axis_2d must be either 'width' or 'length', not {ip_axis_2d}")

    # Note: Often X & Y are used as the distance from the edge to the applied load. X = B/2 - e_B
    # Here e is used to quantify offset vertical load, hload at height and direct moment
    temp_width_eff = temp_fd_width - 2 * abs(e_width)
    temp_length_eff = temp_fd_length - 2 * abs(e_length)
    if temp_length_eff > temp_width_eff:
        length_eff = temp_length_eff
        width_eff = temp_width_eff
        rev = 0
    else:
        length_eff = temp_width_eff
        width_eff = temp_length_eff
        rev = 1

    area_foundation = length_eff * width_eff
    if ip_axis_2d is not None:
        area_foundation = width_eff

    # check CL 3.4.1
    if temp_width_eff / 2 < temp_fd_width / 6 or temp_length_eff / 2 < temp_fd_length / 6:
        print("failed on eccentricity")

    # LOAD FACTORS:
    fd.nq_factor = ((np.tan(np.pi / 4 + sl.phi_r / 2)) ** 2 * np.exp(np.pi * np.tan(sl.phi_r)))
    if sl.phi_r == 0:
        fd.nc_factor = 5.14
    else:
        fd.nc_factor = (fd.nq_factor - 1) / np.tan(sl.phi_r)
    fd.ng_factor = 2.0 * (fd.nq_factor - 1) * np.tan(sl.phi_r)

    # shape factors:
    if length_eff / width_eff > 10:
        s_c = 1.0
        s_q = 1.0
        s_g = 1.0
    else:
        s_c = 1.0 + fd.nq_factor / fd.nc_factor * width_eff / length_eff
        s_q = 1 + width_eff / length_eff * np.tan(sl.phi_r)
        s_g = max(1.0 - 0.4 * width_eff / length_eff, 0.6)  # add limit of 0.6 based on Vesic

    # depth factors:
    if fd.depth / width_eff > 1:
        k = np.arctan(fd.depth / width_eff)
    else:
        k = fd.depth / width_eff
    if sl.phi_r == 0:
        d_c = 1 + 0.4 * k
        d_q = 1.0
    else:
        d_q = (1 + 2 * np.tan(sl.phi_r) * (1 - np.sin(sl.phi_r)) ** 2 * k)
        d_c = d_q - (1 - d_q) / (fd.nq_factor * np.tan(sl.phi_r))
    d_g = 1.0

    # load inclination factors:
    if sl.phi_r == 0:
        i_c = 0.5 * (1 + np.sqrt(1 - horizontal_load / (area_foundation * sl.cohesion)))
        i_q = 1.0
        i_g = 1.0
    else:
        if hload_b == 0:
            i_q = 1 - horizontal_load / (nload + area_foundation * sl.cohesion /
                                         np.tan(sl.phi_r))
            i_g = i_q
        elif hload_b > 0 and hload_l == 0:
            i_q = ((1 - 0.7 * horizontal_load / (nload + area_foundation * sl.cohesion /
                                                 np.tan(sl.phi_r))) ** 3)
            i_g = ((1 - horizontal_load / (nload + area_foundation * sl.cohesion /
                                           np.tan(sl.phi_r))) ** 3)
        else:
            raise DesignError("not setup for bi-directional loading")
        i_c = (i_q * fd.nq_factor - 1) / (fd.nq_factor - 1)

    # ground slope factors:
    g_c = 1 - slope * (1.0 - fd.depth / (2 * width_eff)) / 150
    g_q = (1 - np.tan(slope * (1 - fd.depth / (2 * width_eff)))) ** 2
    g_g = g_q

    # stress at footing base:
    if gwl == 0:
        q_d = sl.unit_bouy_weight * fd.depth
        unit_weight = sl.unit_bouy_weight
    elif 0 < gwl < fd.depth:
        q_d = (sl.unit_dry_weight * gwl) + (sl.unit_bouy_weight * (fd.depth - gwl))
        unit_weight = sl.unit_bouy_weight
    elif fd.depth <= gwl <= fd.depth + width_eff:
        sl.average_unit_bouy_weight = sl.unit_bouy_weight + (
                ((gwl - fd.depth) / width_eff) * (sl.unit_dry_weight - sl.unit_bouy_weight))
        q_d = sl.unit_dry_weight * fd.depth
        unit_weight = sl.average_unit_bouy_weight
    elif gwl > fd.depth + width_eff:
        q_d = sl.unit_dry_weight * fd.depth
        unit_weight = sl.unit_dry_weight
    else:
        raise DesignError(f'gwl must be positive or zero, not: {gwl}')
    q_d += ob  # add overburden

    if verbose:
        log("Nc: %.3f" % fd.nc_factor)
        log("Nq: %.3f" % fd.nq_factor)
        log("Ng: %.3f" % fd.ng_factor)
        log("H: %.3f" % horizontal_load)
        log("s_c: %.3f" % s_c)
        log("s_q: %.3f" % s_q)
        log("s_g: %.3f" % s_g)
        log("d_c: %.3f" % d_c)
        log("d_q: %.3f" % d_q)
        log("d_g: %.3f" % d_g)
        log("i_c: %.3f" % i_c)
        log("i_q: %.3f" % i_q)
        log("i_g: %.3f" % i_g)
        log("g_c: %.3f" % g_c)
        log("g_q: %.3f" % g_q)
        log("g_g: %.3f" % g_g)
        log("q_d: %.3f" % q_d)
    if save_factors:
        fd.width_eff = width_eff
        fd.length_eff = length_eff
        fd.s_c = s_c
        fd.s_q = s_q
        fd.s_g = s_g
        fd.d_c = d_c
        fd.d_q = d_q
        fd.d_g = d_g
        fd.i_c = i_c
        fd.i_q = i_q
        fd.i_g = i_g

        fd.q_d = q_d

    # Capacity
    fd.q_ult = (sl.cohesion * fd.nc_factor * s_c * d_c * i_c * g_c +
                q_d * fd.nq_factor * s_q * d_q * i_q * g_q +
                0.5 * width_eff * unit_weight *
                fd.ng_factor * s_g * d_g * i_g * g_g)
    if rev:
        fd.width_eff = length_eff
        fd.length_eff = width_eff
    else:
        fd.width_eff = width_eff
        fd.length_eff = length_eff
    if ip_axis_2d is not None:
        setattr(fd, f'{fd.ip_axis}_eff', 1)
    fd.area_eff = fd.width_eff * fd.length_eff
    fd.n_ult = fd.q_ult * fd.area_eff
    if verbose:
        log("q_ult: %.3f" % fd.q_ult)
    return fd.q_ult


def capacity_salgado_2008(sl, fd, gwl=1e6, verbose=0, save_factors=0, **kwargs):
    """
    Calculates the capacity according to
     The Engineering of Foundations textbook by Salgado

     ISBN: 0072500581

     The method combines load factors from Bolton (1979) and shape factors from Lyamin et al. (2006)

    Note: if using BH1979 then phi should be from plane strain
    Else: Phi should be appropriate for the shape (not plane strain phi, i.e. should be lower if square than strip)

    :param sl: Soil object
    :param fd: Foundation object
    :param e_length:
        Vertical load eccentricity in length dir
    :param e_width:
        Vertical load eccentricity in width dir

    :param verbose: verbosity
    :return: ultimate bearing stress
    """
    # Need to make adjustments if sand  has DR<40% or
    # clay has liquidity indices greater than 0.7
    if not kwargs.get("disable_requires", False):
        models.check_required(sl, ["phi_r", "cohesion", "unit_dry_weight"])
        models.check_required(fd, ["length", "width", "depth"])
    use_bh1970_factors = kwargs.get('use_bh1970_factors', 0)
    use_loukidis_and_salgado_2011 = kwargs.get('use_loukidis_and_salgado_2011', 0)

    hload_b = kwargs.setdefault('hload_length', 0)
    hload_l = kwargs.setdefault('hload_width', 0)
    nload = kwargs.setdefault('nload', 1)

    # eccentricity due to offset vertical load, hload at height and direct moment.
    e_length = kwargs.setdefault('e_length', 0)
    e_width = kwargs.setdefault('e_width', 0)

    ip_axis_2d = kwargs.get('ip_axis_2d', None)
    ob = kwargs.get('ob', 0.0)  # overburden pressure

    temp_fd_length = fd.length
    temp_fd_width = fd.width
    if ip_axis_2d is not None:
        if ip_axis_2d == 'width':
            temp_fd_length = temp_fd_width * 100
            e_length = temp_fd_length / 2
        elif ip_axis_2d == 'length':
            temp_fd_width = temp_fd_length * 100
            e_width = temp_fd_width / 2
        else:
            raise ValueError(f"ip_axis_2d must be either 'width' or 'length', not {ip_axis_2d}")

    temp_width_eff = temp_fd_width - 2 * abs(e_width)
    temp_length_eff = temp_fd_length - 2 * abs(e_length)
    if temp_length_eff > temp_width_eff:
        length_eff = temp_length_eff
        width_eff = temp_width_eff
        rev = 0
    else:
        length_eff = temp_width_eff
        width_eff = temp_length_eff
        rev = 1

    # LOAD FACTORS:
    fd.nq_factor = np.exp(np.pi * np.tan(sl.phi_r)) * (1 + np.sin(sl.phi_r)) / (1 - np.sin(sl.phi_r))  # Eq 10.6
    fd.ng_factor = (fd.nq_factor - 1) * np.tan(1.32 * sl.phi_r)  # Eq 10.13
    if use_bh1970_factors:
        fd.ng_factor = 1.5 * (fd.nq_factor - 1) * np.tan(sl.phi_r)  # BH1970 (Eq 10.12)

    if use_loukidis_and_salgado_2011:
        fd.ng_factor = (fd.nq_factor - 0.6) * np.tan(1.33 * sl.phi_r)
    if sl.phi_r == 0:
        fd.nc_factor = 5.14
    else:
        fd.nc_factor = (fd.nq_factor - 1) / np.tan(sl.phi_r)  # Eq 10.11 (application of Coquot's principle)
    d_o_b_min = 2.0  # Range from Lyamin is 0-2
    # shape factors:
    if ip_axis_2d is not None:
        s_q = 1.0
        s_g = 1.0
        s_c = 1.0
    else:
        if use_bh1970_factors:
            s_q = 1 + (width_eff / length_eff) * np.sin(sl.phi_r)
            s_g = max(1 - 0.4 * width_eff / length_eff, 0.6)
            if sl.phi_r == 0:
                s_c = 1 + 0.2 * width_eff / length_eff
            else:
                s_c = 1.0 + fd.nq_factor / fd.nc_factor * width_eff / length_eff
        else:
            # From Table 11-8 on page 496, note phi is in degrees
            d_o_b = min(fd.depth / width_eff, d_o_b_min)
            b_o_l = width_eff / length_eff
            if d_o_b == 0:
                s_q = 1.0  # note that equation does tend to 1 as d_o_b -> 0
            else:
                s_q = 1 + (0.0952 * sl.phi - 1.60) * d_o_b ** (0.583 - 0.0079 * sl.phi) * b_o_l ** (1 - 0.15 * d_o_b)
            s_g = 1 + (0.0345 * sl.phi -1.0611) * b_o_l  # note: in original paper Lyamin 2007, Eq is 1 + (0.0336*phi - 1) * b_o_l
            s_c = calc_shape_factor_coh_salgado_et_al_2004(fd)  # Salgado (2008) Eq. 10.22

    # depth factors:
    d_o_b = min(d_o_b_min, fd.depth / width_eff)

    if use_bh1970_factors:
        d_q = 1 + 2 * np.tan(sl.phi_r) * (1 - np.sin(sl.phi_r)) ** 2 * d_o_b
    else:
        d_o_bm = 0.01  # no limit set in paper but fits well to FEM results and matches well to figure 11 in Lyamin
        if d_o_b <= d_o_bm:
            # d_q = 3.0  # note that it is inverse to d/b see paper Lyamin et al. (2007) Eq 14
            d_q = 1 + (0.0044 * sl.phi + 0.356) * d_o_bm ** -0.28
        else:
            d_q = 1 + (0.0044 * sl.phi + 0.356) * d_o_b ** -0.28  # Lyamin et al. (2006) [Salgado Table 10-7]
    d_g = 1.0  # Both Lyamin and BH1970 (Table 10-7)
    d_c = 1.0 + 0.27 * np.sqrt(d_o_b)

    # stress at footing base:
    if gwl == 0:
        q_d = sl.unit_bouy_weight * fd.depth
        unit_weight = sl.unit_bouy_weight
    elif 0 < gwl < fd.depth:
        q_d = (sl.unit_dry_weight * gwl) + (sl.unit_bouy_weight * (fd.depth - gwl))
        unit_weight = sl.unit_bouy_weight
    elif fd.depth <= gwl <= fd.depth + width_eff:
        sl.average_unit_bouy_weight = sl.unit_bouy_weight + (
                ((gwl - fd.depth) / width_eff) * (sl.unit_dry_weight - sl.unit_bouy_weight))
        q_d = sl.unit_dry_weight * fd.depth
        unit_weight = sl.average_unit_bouy_weight
    elif gwl > fd.depth + width_eff:
        q_d = sl.unit_dry_weight * fd.depth
        unit_weight = sl.unit_dry_weight
    else:
        raise DesignError(f'gwl must be positive or zero, not: {gwl}')
    q_d += ob  # add overburden

    if verbose:
        log("width_eff: ", width_eff)
        log("length_eff: ", length_eff)
        log("Nc: ", fd.nc_factor)
        log("Nq: ", fd.nq_factor)
        log("Ng: ", fd.ng_factor)
        log("s_c: %.3f" % s_c)
        log("s_q: %.3f" % s_q)
        log("s_g: %.3f" % s_g)
        log("d_c: %.3f" % d_c)
        log("d_q: %.3f" % d_q)
        log("d_g: %.3f" % d_g)
        log("q_d: %.3f" % q_d)
    if save_factors:
        fd.width_eff = width_eff
        fd.length_eff = length_eff
        fd.s_c = s_c
        fd.s_q = s_q
        fd.s_g = s_g
        fd.d_c = d_c
        fd.d_q = d_q
        fd.d_g = d_g
        fd.q_d = q_d

    # Capacity
    fd.q_ult = (sl.cohesion * fd.nc_factor * s_c * d_c +
                q_d * fd.nq_factor * s_q * d_q +
                0.5 * width_eff * unit_weight *
                fd.ng_factor * s_g * d_g)

    if rev:
        fd.width_eff = length_eff
        fd.length_eff = width_eff
    else:
        fd.width_eff = width_eff
        fd.length_eff = length_eff
    if ip_axis_2d is not None:
        setattr(fd, f'{fd.ip_axis}_eff', 1)
    fd.area_eff = fd.width_eff * fd.length_eff
    fd.n_ult = fd.q_ult * fd.area_eff
    if verbose:
        log("q_ult: %.3f" % fd.q_ult)
    return fd.q_ult


def size_footing_for_capacity(sl, vertical_load, fos=1.0, length_to_width=1.0, verbose=0, unit_weight=0, **kwargs):
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
        fs_actual = bearing_capacity / (vertical_load + unit_weight * fd.area * fd.depth)

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

        fs_array.append(capacity / (vertical_load + unit_weight * fd.area * fd.depth))

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


def calc_crit_span(sl, fd, vertical_load, ip_axis=None, verbose=0, ip_axis_2d=None, method='vesic', ob=0.0):
    """
    Determine the size of a footing given an aspect ratio and a load

    :param sl: Soil object
    :param vertical_load: The applied load to the foundation
    :param fos: The target factor of safety
    :param length_to_width: The desired length to width ratio of the foundation
    :param verbose: verbosity
    :return: a Foundation object
    """

    # Find approximate size
    new_fd = models.RaftFoundation()
    new_fd.width = fd.width
    new_fd.depth = fd.depth
    new_fd.length = fd.length
    if ip_axis is None:
        ip_axis = fd.ip_axis
    prev_ub_len = getattr(fd, ip_axis)
    q_ult = capacity_method_selector(sl, new_fd, method, verbose=max(0, verbose - 1), ip_axis_2d=ip_axis_2d, ob=ob)
    if ip_axis_2d is None:
        area = fd.area
    else:
        area = getattr(fd, ip_axis_2d)
    init_fos = (q_ult * area) / vertical_load
    if init_fos < 1.0:
        raise ValueError
    prev_lb_len = 0  # should be FOS lower than 1.
    l_ip = getattr(fd, ip_axis)
    est_len = l_ip / init_fos
    prev_q = q_ult
    for i in range(50):
        setattr(new_fd, ip_axis, est_len)
        q = capacity_method_selector(sl, new_fd, method, verbose=max(0, verbose - 1), ip_axis_2d=ip_axis_2d, ob=ob)
        if ip_axis_2d is None:
            area = new_fd.area
        else:
            area = getattr(new_fd, ip_axis_2d)
        curr_fos = (q * area) / vertical_load
        if np.isclose(curr_fos, 1.0, rtol=0.01):
            break
        elif curr_fos < 1:
            prev_lb_len = est_len
            est_len = (prev_ub_len + est_len) / 2
        else:
            prev_ub_len = est_len
            est_len = (prev_lb_len + est_len) / 2
    if i == 99:
        raise ValueError(init_fos, curr_fos, est_len, prev_lb_len, prev_ub_len)
    return est_len


def calc_crit_span_and_phi_p_via_salgado_2008(sl, fd, vertical_load, ip_axis=None, verbose=0, ip_axis_2d=None,
                                              **kwargs):
    """
    Determine the size of a footing given an aspect ratio and a load

    :param sl: Soil object
    :param vertical_load: The applied load to the foundation
    :param fos: The target factor of safety
    :param length_to_width: The desired length to width ratio of the foundation
    :param verbose: verbosity
    :return: a Foundation object
    """
    sl_org_phi = sl.phi
    # Find approximate size
    new_fd = sm.RaftFoundation()
    new_fd.width = fd.width
    new_fd.depth = fd.depth
    new_fd.length = fd.length
    if ip_axis_2d:
        ip_axis = ip_axis_2d
    elif ip_axis is None:
        ip_axis = fd.ip_axis
    prev_ub_len = getattr(fd, ip_axis)
    p_eff = calc_m_eff_bc_via_loukidis_and_salgado_2006(sl, fd, ip_axis_2d=ip_axis_2d)
    if ip_axis_2d is None:
        l_o_b = new_fd.llong / new_fd.lshort
    else:
        l_o_b = 7
    pkwargs = {}
    if hasattr(sl, 'phi_c_ps_diff'):
        pkwargs['phi_c_ps_diff'] = sl.phi_c_ps_diff
    sl.phi = calc_phi_peak_fd_salgado_2008(sl.phi_c_txc, p_eff, sl.relative_density, l_o_b, **pkwargs)

    q_ult = capacity_salgado_2008(sl, new_fd, verbose=max(0, verbose - 1), ip_axis_2d=ip_axis_2d, **kwargs)
    if ip_axis_2d is None:
        area = fd.area
    else:
        area = getattr(fd, ip_axis_2d)
    init_fos = (q_ult * area) / vertical_load
    if init_fos < 1.0:
        raise ValueError
    prev_lb_len = 0  # should be FOS lower than 1.
    l_ip = getattr(fd, ip_axis)
    est_len = l_ip / init_fos
    prev_q = q_ult
    for i in range(50):
        setattr(new_fd, ip_axis, est_len)

        p_eff = calc_m_eff_bc_via_loukidis_and_salgado_2006(sl, new_fd, ip_axis_2d=ip_axis_2d, **kwargs)
        if ip_axis_2d is None:
            l_o_b = new_fd.llong / new_fd.lshort
        else:
            l_o_b = 7
        pkwawrgs = {}
        if hasattr(sl, 'phi_c_ps_diff'):
            pkwawrgs['phi_c_ps_diff'] = sl.phi_c_ps_diff
        sl.phi = calc_phi_peak_fd_salgado_2008(sl.phi_c_txc, p_eff, sl.relative_density, l_o_b, **pkwawrgs)
        q = capacity_salgado_2008(sl, new_fd, verbose=max(0, verbose - 1), ip_axis_2d=ip_axis_2d, **kwargs)
        if ip_axis_2d is None:
            area = new_fd.area
        else:
            area = getattr(new_fd, ip_axis_2d)
        curr_fos = (q * area) / vertical_load
        if np.isclose(curr_fos, 1.0, rtol=0.01):
            break
        elif curr_fos < 1:
            prev_lb_len = est_len
            est_len = (prev_ub_len + est_len) / 2
        else:
            prev_ub_len = est_len
            est_len = (prev_lb_len + est_len) / 2
    if i == 99:
        raise ValueError(init_fos, curr_fos, est_len, prev_lb_len, prev_ub_len)
    phi_p = sl.phi
    sl.phi = sl_org_phi
    return est_len, phi_p


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
    elif method == 'nzs':
        return capacity_nzs_vm4_2011(sl, fd, **kwargs)
    elif method == 'terzaghi':
        return capacity_terzaghi_1943(sl, fd, **kwargs)
    elif method == 'brinch_hansen':
        return capacity_brinch_hansen_1970(sl, fd, **kwargs)
    elif method == 'meyerhof':
        return capacity_meyerhof_1963(sl, fd, **kwargs)
    elif method == 'salgado':
        return capacity_salgado_2008(sl, fd, **kwargs)
    else:
        raise ValueError(f"{method} not found. method must be 'vesic', 'nzs', 'terzaghi', 'meyerhof', 'salgado', 'brinch_hansen'")


available_methods = {
    "nzs": capacity_nzs_vm4_2011,
    "terzaghi": capacity_terzaghi_1943,
    "brinch_hansen": capacity_brinch_hansen_1970,
    "meyerhof": capacity_meyerhof_1963,
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


def capacity_meyerhof_and_hanna_1978(sl_0, sl_1, h0, fd, gwl=1e6, ip_axis_2d=None, verbose=0):
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
    return capacity_sp_meyerhof_and_hanna_1978(sp, fd, ip_axis_2d=ip_axis_2d)


def capacity_sp_meyerhof_and_hanna_1978(sp, fd, ip_axis_2d=None, verbose=0):
    """
    Calculates the two-layered foundation capacity according Meyerhof and Hanna (1978)

    :param sp: Soil profile object
    :param fd: Foundation object
    :param wtl: water table level
    :param verbose: verbosity
    :return: ultimate bearing stress
    """
    assert isinstance(sp, sm.SoilProfile)
    temp_fd_length = fd.length
    temp_fd_width = fd.width
    if ip_axis_2d is not None:
        if ip_axis_2d == 'width':
            temp_fd_length = temp_fd_width * 100
        elif ip_axis_2d == 'length':
            temp_fd_width = temp_fd_length * 100
        else:
            raise ValueError(f"ip_axis_2d must be either 'width' or 'length', not {ip_axis_2d}")
    if temp_fd_length > temp_fd_width:  # TODO: deal with plane strain
        fd_length = temp_fd_length
        fd_width = temp_fd_width
    else:
        fd_length = temp_fd_width
        fd_width = temp_fd_length

    sl_1 = sp.layer(1)
    sl_2 = sp.layer(2)
    h0 = sp.layer_depth(2)
    gwl = sp.gwl

    sl_1.nq_factor_1 = (
                (np.tan(np.pi / 4 + np.deg2rad(sl_1.phi / 2))) ** 2 * np.exp(np.pi * np.tan(np.deg2rad(sl_1.phi))))
    if sl_1.phi == 0:
        sl_1.nc_factor_1 = 5.14
    else:
        sl_1.nc_factor_1 = (sl_1.nq_factor_1 - 1) / np.tan(np.deg2rad(sl_1.phi))
    sl_1.ng_factor_1 = (sl_1.nq_factor_1 - 1) * np.tan(1.4 * np.deg2rad(sl_1.phi))

    sl_2.nq_factor_2 = (
            (np.tan(np.pi / 4 + np.deg2rad(sl_2.phi / 2))) ** 2 * np.exp(np.pi * np.tan(np.deg2rad(sl_2.phi))))
    if sl_2.phi == 0:
        sl_2.nc_factor_2 = 5.14
    else:
        sl_2.nc_factor_2 = (sl_2.nq_factor_2 - 1) / np.tan(np.deg2rad(sl_2.phi))
    sl_2.ng_factor_2 = (sl_2.nq_factor_2 - 1) * np.tan(1.4 * np.deg2rad(sl_2.phi))

    if verbose:
        log("Nc: ", sl_2.nc_factor_2)
        log("Nq: ", sl_2.nq_factor_2)
        log("Ng: ", sl_2.ng_factor_2)

    sl_1.kp_1 = (np.tan(np.pi / 4 + np.deg2rad(sl_1.phi / 2))) ** 2
    sl_2.kp_2 = (np.tan(np.pi / 4 + np.deg2rad(sl_2.phi / 2))) ** 2

    # shape factors
    if sl_1.phi >= 10:
        sl_1.s_c_1 = 1 + 0.2 * sl_1.kp_1 * (fd_width / fd_length)
        sl_1.s_q_1 = 1.0 + 0.1 * sl_1.kp_1 * (fd_width / fd_length)
    else:
        sl_1.s_c_1 = 1 + 0.2 * (fd_width / fd_length)
        sl_1.s_q_1 = 1.0
    sl_1.s_g_1 = sl_1.s_q_1

    if sl_2.phi >= 10:
        sl_2.s_c_2 = 1 + 0.2 * sl_2.kp_2 * (fd_width / fd_length)
        sl_2.s_q_2 = 1.0 + 0.1 * sl_2.kp_2 * (fd_width / fd_length)
    else:
        sl_2.s_c_2 = 1 + 0.2 * (fd_width / fd_length)
        sl_2.s_q_2 = 1.0
    sl_2.s_g_2 = sl_2.s_q_2

    # Note: this method explicitly accounts for the foundation depth, so there are no depth factors
    # TODO: inclination factors, see doi.org/10.1139/t78-060

    # Capacity
    s = 1

    if ip_axis_2d is not None:
        r = 1
    else:
        r = 1 + (fd_width / fd_length)

    # put the same things before that condition
    # effective weight not in the soil object

    if gwl == 0:  # case 1: GWL at surface
        q_at_interface = sl_1.unit_bouy_weight * h0
        unit_eff_weight_1_at_fd_depth = sl_1.unit_bouy_weight
        unit_eff_weight_1_at_interface = sl_1.unit_bouy_weight
        unit_eff_weight_2_below_foundation = sl_2.unit_bouy_weight

    elif 0 < gwl <= fd.depth:  # Case 2: GWL at between foundation depth and surface
        q_at_interface = (sl_1.unit_dry_weight * gwl) + (sl_1.unit_bouy_weight * (h0 - gwl))
        q_d = (sl_1.unit_dry_weight * gwl) + (sl_1.unit_bouy_weight * (fd.depth - gwl))
        unit_eff_weight_1_at_fd_depth = q_d / fd.depth
        unit_eff_weight_1_at_interface = sl_1.unit_bouy_weight
        unit_eff_weight_2_below_foundation = sl_2.unit_bouy_weight

    elif fd.depth < gwl <= fd_width + fd.depth:
        if gwl < h0:  # Case 3: GWL at between foundation depth and foundation depth plus width, and GWL < layer 1 depth

            average_unit_bouy_weight = sl_1.unit_bouy_weight + (
                    ((gwl - fd.depth) / fd_width) * (sl_1.unit_dry_weight - sl_1.unit_bouy_weight))

            q_at_interface = (sl_1.unit_dry_weight * gwl) + (sl_1.unit_bouy_weight * (h0 - gwl))
            unit_eff_weight_1_at_fd_depth = sl_1.unit_dry_weight
            unit_eff_weight_1_at_interface = average_unit_bouy_weight
            unit_eff_weight_2_below_foundation = sl_2.unit_bouy_weight

        else:  # Case 4: GWL at between foundation depth and foundation depth plus width, and GWL > layer 1 depth
            average_unit_bouy_weight = sl_2.unit_bouy_weight + (
                    ((gwl - h0) / fd_width) * (sl_2.unit_dry_weight - sl_2.unit_bouy_weight))

            q_at_interface = sl_1.unit_dry_weight * h0
            unit_eff_weight_1_at_fd_depth = sl_1.unit_dry_weight
            unit_eff_weight_1_at_interface = sl_1.unit_dry_weight
            unit_eff_weight_2_below_foundation = average_unit_bouy_weight

    elif gwl > fd.depth + fd_width:  # Case 5: GWL beyond foundation depth plus width
        q_at_interface = sl_1.unit_dry_weight * h0
        unit_eff_weight_1_at_fd_depth = sl_1.unit_dry_weight
        unit_eff_weight_1_at_interface = sl_1.unit_dry_weight
        unit_eff_weight_2_below_foundation = sl_2.unit_dry_weight
    else:
        raise ValueError("Could not interpret inputs")  # never reached

    # Das Eq 4.33 (q1)
    q_1 = (sl_1.cohesion * sl_1.nc_factor_1) + (0.5 * unit_eff_weight_1_at_interface * fd_width * sl_1.ng_factor_1)
    # Das Eq 4.34 (q2)
    q_2 = (sl_2.cohesion * sl_2.nc_factor_2) + (0.5 * unit_eff_weight_2_below_foundation * fd_width * sl_2.ng_factor_2)

    q_ult5 = r * (unit_eff_weight_1_at_interface * ((h0 - fd.depth) ** 2)) * (1 + (2 * fd.depth / (h0 - fd.depth))) * (
                np.tan(np.deg2rad(sl_1.phi)) / fd_width) * s

    q1_q0 = q_2 / q_1

    # calculate the ca factor
    # if sl_1.cohesion == 0:
    #     c1_c0 = 0
    # else:
    #     c1_c0 = sl_2.cohesion / sl_1.cohesion
    x = np.array([0.000, 0.082, 0.206, 0.298, 0.404, 0.509, 0.598, 0.685, 0.772])
    y = np.array([0.627, 0.700, 0.794, 0.855, 0.912, 0.948, 0.968, 0.983, 0.997])

    # raise Warning("ca should be interpolated using q1/q2 not cohesion, see Figure 4 in MH1978")
    ca_c0 = np.interp(q1_q0, x, y)
    ca = ca_c0 * sl_1.cohesion  # Eq 2?

    # ks - coefficient of punching shear
    x_0 = np.array([0, 20.08, 22.42, 25.08, 27.58, 30.08, 32.58, 34.92, 37.83, 40.00, 42.67, 45.00, 47.00, 49.75])
    y_0 = np.array([0.93, 0.93, 0.93, 0.93, 1.01, 1.17, 1.32, 1.56, 1.87, 2.26, 2.72, 3.35, 3.81, 4.82])
    x_2 = np.array([0, 20.08, 22.50, 25.08, 27.58, 30.08, 32.50, 35.00, 37.67, 40.17, 42.67, 45.00, 47.50, 50.00])
    y_2 = np.array([1.55, 1.55, 1.71, 1.86, 2.10, 2.33, 2.72, 3.11, 3.81, 4.43, 5.28, 6.14, 7.46, 9.24])
    x_4 = np.array([0, 20.00, 22.51, 25.10, 27.69, 30.11, 32.45, 35.04, 37.88, 40.14, 42.65, 45.07, 47.33, 50.08])
    y_4 = np.array([2.49, 2.49, 2.64, 2.87, 3.34, 3.81, 4.43, 5.20, 6.29, 7.38, 9.01, 11.11, 14.29, 19.34])
    x_10 = np.array([0, 20.00, 22.50, 25.08, 28.00, 30.00, 32.50, 34.92, 37.50, 40.17, 42.42, 45.00, 47.17, 50.08])
    y_10 = np.array([3.27, 3.27, 3.74, 4.44, 5.37, 6.07, 7.16, 8.33, 10.04, 12.30, 15.95, 21.17, 27.47, 40.00])
    x_int = sl_1.phi

    if sl_1.phi < 1:
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
            fd.ks = None

            # raise DesignError(
            #     "Cannot compute 'ks', bearing ratio out-of-range (q1_q0 = %.3f) required: 0-1." % q1_q0)

    # qb  # MH78 Eq 6
    q_b1 = (sl_2.cohesion * sl_2.nc_factor_2 * sl_2.s_c_2)
    q_b2 = (q_at_interface * sl_2.nq_factor_2 * sl_2.s_q_2)
    q_b3 = (unit_eff_weight_2_below_foundation * fd_width * sl_2.ng_factor_2 * sl_2.s_g_2 / 2)
    q_b = q_b1 + q_b2 + q_b3

    # Das Eq 4.35 or MH78 Eq 7
    q_t1 = (sl_1.cohesion * sl_1.nc_factor_1 * sl_1.s_c_1)
    q_t2 = (unit_eff_weight_1_at_fd_depth * fd.depth * sl_1.nq_factor_1 * sl_1.s_q_1)
    q_t3 = (unit_eff_weight_1_at_interface * fd_width * sl_1.ng_factor_1 * sl_1.s_g_1 / 2)
    q_t = q_t1 + q_t2 + q_t3

    if fd.ks is not None:  # bottom soil is weaker
        # qu  # MH78 Eq 5
        a = 1  # assumed to  be one but can range between 1.1 and 1.27 for square footings according to Das (1999) Ch 4
        q_ult4 = (r * (2 * ca * (h0 - fd.depth) / fd_width) * a)  # Note: slightly different to MH78
        q_ult5_ks = q_ult5 * fd.ks
        q_ult6 = q_at_interface - unit_eff_weight_1_at_fd_depth * fd.depth
        q_ult = q_b + q_ult4 + q_ult5_ks - q_ult6  # Eq 5 from MH78 (Das Eq 4.36)

        if q_ult > q_t:
            if h0 > fd_width / 2:
                fd.q_ult = q_t

            else:
                vert_eff_stress_interface = sp.get_v_eff_stress_at_depth(h0)
                vert_eff_stress_lowest = sp.get_v_eff_stress_at_depth(fd_width + fd.depth)
                average_eff_stress = (vert_eff_stress_interface + vert_eff_stress_lowest) / 2

                c_2_eff = sl_2.cohesion + average_eff_stress * np.tan(np.radians(sl_2.phi))

                if sl_1.cohesion > c_2_eff:
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
                    c1_over_c2 = sl_1.cohesion / c_2_eff

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
                        m = np.interp(h_over_b, [0.1, 0.125], [m_a, m_b])
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

                    fd.q_ult = (sl_1.cohesion * m * sl_1.nc_factor_1) + (unit_eff_weight_1_at_fd_depth * fd.depth)
        else:
            fd.q_ult = q_ult
    else:  # top soil is weaker
        if sl_1.cohesion > 0:
            h_f = 1.0 * fd_width
        else:
            phi_loose = 30.0
            phi_dense = 42.0
            h_f = np.clip(1 + (sl_1.phi - phi_loose) / (phi_dense - phi_loose), 1, 2) * fd_width
        fd.q_ult = max([q_t + (q_b - q_t) * max([(1 - (h0 - fd.depth) / h_f), 0]) ** 2, q_t])

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