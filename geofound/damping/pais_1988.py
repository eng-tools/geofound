import numpy as np

import geofound as gf


def calc_rot_via_pais_1988(sl, fd, ip_axis=None, a0=0.0, **kwargs):
    """
    Rotational damping of foundation from Pais and Kausel (1988).

    Parameters
    ----------
    sl: Soil Object
    fd: Foundation object
    ip_axis: str
        The axis that is in the plane of deformation
    a0: float
        dynamic factor

    Returns
    -------
    beta_f: float
        Rotational damping of the foundation
    """
    if not kwargs.get("disable_requires", False):
        gf.models.check_required(sl, ["g_mod", "poissons_ratio"])
        gf.models.check_required(fd, ["length", "width", "depth"])

    if fd.i_ww >= fd.i_ll:
        len_dominant = True
        l = fd.length * 0.5
        b = fd.width * 0.5
    else:
        len_dominant = False
        l = fd.width * 0.5
        b = fd.length * 0.5
    if ip_axis is None:
        ip_axis = fd.ip_axis

    if (ip_axis == 'width' and len_dominant) or (ip_axis == 'length' and not len_dominant):
        xx_axis = True  # weaker rotation (rotation about longest length)
    else:
        xx_axis = False

    v = sl.poissons_ratio
    chi = min(np.sqrt(2 * (1 - v) / (1 - 2 * v)), 2.5)
    l_o_b = l / b
    if xx_axis:
        k_xx_static = gf.stiffness.calc_rot_via_pais_1988(sl, fd, ip_axis=ip_axis, a0=0)
        k_xx_dyn = gf.stiffness.calc_rot_via_pais_1988(sl, fd, ip_axis=ip_axis, a0=a0)
        alp_xx = k_xx_dyn / k_xx_static
        if fd.depth is not None and fd.depth != 0.0:
            d = fd.depth
            d_o_b = d / b
            if fd.depth < 0.0:
                raise ValueError(f'foundation depth must be zero or greater, not {fd.depth}')
            # Table 2-3b (NIST 2013)
            beta = ((4. / 3) * (d_o_b + d_o_b ** 3 + chi * l_o_b * d_o_b ** 3 + 3 * d_o_b * l_o_b + chi * l_o_b) * a0 ** 2 /
                    (k_xx_static / (sl.g_mod * b ** 3)) * ((1.8 / (1 + 1.75 * (l_o_b - 1)) + a0 ** 2)) +
                    (4. / 3) * (chi * l / b + 1) * d_o_b ** 3 /
                    (k_xx_static / (sl.g_mod * b ** 3))) * (
                              a0 / (2 * alp_xx))
        else:
            # Table 2-3a
            beta = ((4 * chi / 3) * l_o_b * a0 ** 2 /
                    ((k_xx_static / (sl.g_mod * b ** 3)) * ((2.2 - 0.4 / l_o_b ** 3) + a0 ** 2))) * \
                   (a0 / 2 / alp_xx)
    else:
        k_yy_static = gf.stiffness.calc_rot_via_pais_1988(sl, fd, ip_axis=ip_axis, a0=0)
        k_yy_dyn = gf.stiffness.calc_rot_via_pais_1988(sl, fd, ip_axis=ip_axis, a0=a0)
        alp_yy = k_yy_dyn / k_yy_static
        # Table 2-3b (NIST 2013)
        if fd.depth is not None and fd.depth != 0.0:
            d = fd.depth
            if fd.depth < 0.0:
                raise ValueError(f'foundation depth must be zero or greater, not {fd.depth}')
            beta = ((4. / 3) * ((l / b) ** 3 * (d / b) + chi * (d / b) ** 3 * (l / b) + (d / b) ** 3 +
                                3 * (d / b) * (l / b) ** 2 + chi * (l / b) ** 3) * a0 ** 2 /
                    (k_yy_static / (sl.g_mod * b ** 3)) * (1.8 / (1 + 1.75 * (l / b - 1)) + a0 ** 2) +
                    (4. / 3) * (l / b + chi) * (d / b) ** 3 /
                    (k_yy_static / (sl.g_mod * b ** 3))) * \
                   (a0 / (2 * alp_yy))
        else:
            beta = ((4 * chi / 3) * (l / b) * a0 ** 2 /
                    k_yy_static / ((sl.g_mod * b ** 3) * (1.8 / (1 + 1.75 * (l_o_b - 1))) + a0 ** 2)) * \
                   (a0 / 2 / alp_yy)

    return beta



def calc_vert_via_pais_1988(sl, fd, a0=0):
    """
    Vertical damping of foundation from Pais and Kausel (1988).

    Parameters
    ----------
    sl: sm.Soil object
        Soil object
    fd: sm.Foundation object
        Foundation object
    a0: float
        Normalised frequency

    Returns
    -------

    """
    v = sl.poissons_ratio
    l = fd.length * 0.5
    b = fd.width * 0.5

    k_n_static = gf.stiffness.calc_vert_via_pais_1988(sl, fd, a0=0)
    k_n_dyn = gf.stiffness.calc_vert_via_pais_1988(sl, fd, a0=a0)
    alp_z = k_n_dyn / k_n_static
    chi = min(np.sqrt(2 * (1 - v) / (1 - 2 * v)), 2.5)
    l_o_b = l / b
    if fd.depth is not None and fd.depth != 0.0:
        if fd.depth < 0.0:
            raise ValueError(f'foundation depth must be zero or greater, not {fd.depth}')
        d_o_b = fd.depth / b
        beta = (4 * (chi * l_o_b + d_o_b * (1 + l_o_b)) / (k_n_static / sl.g_mod / b)) * (a0 / (2 * alp_z))
    else:
        beta = (4 * chi * l_o_b / (k_n_static / sl.g_mod / b)) * (a0 / (2 * alp_z))
    return beta


def calc_horz_via_pais_1988(sl, fd, ip_axis=None, a0=0.0):
    """
    Calculate the horizontal damping for translation along an axis from Pais and Kausel (1988).

    Parameters
    ----------
    sl: Soil Object
    fd: Foundation object
    ip_axis: str
        The axis that is in the plane of deformation
    a0: float
        dynamic factor

    Returns
    -------

    """
    if fd.length >= fd.width:
        len_dominant = True
        l = fd.length * 0.5
        b = fd.width * 0.5
    else:
        len_dominant = False
        l = fd.width * 0.5
        b = fd.length * 0.5
    if ip_axis is None:
        ip_axis = fd.ip_axis
    if ip_axis is None:
        raise ValueError('ip_axis must be set')

    if (ip_axis == 'length' and len_dominant) or (ip_axis == 'width' and not len_dominant):
        x_axis = True  # Direction of l
    else:
        x_axis = False  # Direction of b
    v = sl.poissons_ratio
    chi = min(np.sqrt(2 * (1 - v) / (1 - 2 * v)), 2.5)
    l_o_b = l / b
    k_h = gf.stiffness.calc_horz_via_pais_1988(sl, fd, ip_axis=ip_axis, a0=a0)
    alp = 1.0  # note alpha_x = alpha_y = 1
    if x_axis:  # x_axis
        if fd.depth:
            d_o_b = fd.depth / b
            beta_h = (4 * (l_o_b + d_o_b * (chi + l_o_b)) /
                      (k_h / sl.g_mod / b)) * \
                     (a0 / (2 * alp))
        else:
            beta_h = ((4 * l_o_b) / (k_h / sl.g_mod / b)) * (a0 / (2 * alp))

    else:
        if fd.depth:
            d_o_b = fd.depth / b
            beta_h = (4 * (l_o_b + d_o_b * (1 + chi * l_o_b)) /
                      (k_h / sl.g_mod / b)) * \
                     (a0 / (2 * alp))
        else:
            beta_h = ((4 * l_o_b) / (k_h / sl.g_mod / b)) * (a0 / (2 * alp))

    return beta_h
