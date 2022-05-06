import numpy as np

import geofound as gf


def calc_rot_via_pais_1988(sl, fd, ip_axis=None, a0=0.0, **kwargs):
    """
    Rotational stiffness of foundation from Pais and Kausel (1988).

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
    k_f: float
        Rotational stiffness of the foundation
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
    if (ip_axis == 'width' and len_dominant) or (ip_axis == 'length' and not len_dominant):
        xx_axis = True  # weaker rotation (rotation about longest length)
    else:
        xx_axis = False

    v = sl.poissons_ratio
    chi = min(np.sqrt(2 * (1 - v) / (1 - 2 * v)), 2.5)
    if xx_axis:
        # x-axis
        # Table 2-3a (NIST 2013)
        f_r_dyn_surf = 1 - ((0.55 + 0.01 * (l / b - 1) ** 0.5) * a0 ** 2 / ((2.4 - 0.4 / (l / b) ** 3) + a0 ** 2))
        k_f_0_static_surf = (sl.g_mod * b ** 3 / (1 - v) * (3.2 * (l / b) + 0.8))
        if fd.depth is not None and fd.depth != 0.0:
            d = fd.depth
            if fd.depth < 0.0:
                raise ValueError(f'foundation depth must be zero or greater, not {fd.depth}')
            n_emb = 1.0 + fd.depth / b + (1.6 / (0.35 + (l / b)) * (fd.depth / b) ** 2)
            k_f_0_static = k_f_0_static_surf * n_emb

        else:
            n_emb = 1
            k_f_0_static = k_f_0_static_surf
    else:
        # Table 2-3a (NIST 2013)
        f_r_dyn_surf = 1 - (0.55 * a0 ** 2 / ((0.6 + 1.4 / (l / b) ** 3) + a0 ** 2))

        k_f_0_static_surf = (sl.g_mod * b ** 3 / (1 - v) * (3.73 * (l / b) ** 2.4 + 0.27))
        if fd.depth is not None and fd.depth != 0.0:
            d = fd.depth
            if fd.depth < 0.0:
                raise ValueError(f'foundation depth must be zero or greater, not {fd.depth}')
            n_emb = 1.0 + d / b + (1.6 / (0.35 + (l / b) ** 4) * (d / b) ** 2)
            k_f_0_static = k_f_0_static_surf * n_emb

        else:
            n_emb = 1
            k_f_0_static = k_f_0_static_surf
    f_r_dyn = f_r_dyn_surf  # Note at Table 2-3b f_dyn_emb = f_dyn_surf
    return k_f_0_static * f_r_dyn * n_emb


def calc_vert_via_pais_1988(sl, fd, a0=0):
    """
    Vertical stiffness of foundation from Pais and Kausel (1988).

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
    k_v_0 = sl.g_mod * b / (1 - v) * (3.1 * (l / b) ** 0.75 + 1.6)
    if a0:  # Table 2-3a
        f_dyn_surf = 1 - ((0.4 + 0.2 / (l / b)) * a0 ** 2 / (10 / (1 + 3 * (l / b) - 1) + a0 ** 2))
    else:
        f_dyn_surf = 1
    if fd.depth is not None and fd.depth != 0.0:
        if fd.depth < 0.0:
            raise ValueError(f'foundation depth must be zero or greater, not {fd.depth}')
        n_emb = 1 + (0.25 + 0.25 / (l / b)) * (fd.depth / b) ** 0.8
    else:
        n_emb = 1
    f_dyn = f_dyn_surf  # See note end of Table 2-3b
    return k_v_0 * n_emb * f_dyn


def calc_horz_via_pais_1988(sl, fd, ip_axis=None, a0=0.0):
    """
    Calculate the horizontal stiffness for translation along an axis.

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
    if ip_axis not in ['length', 'width']:
        raise ValueError(f'Must set ip_axis to either "length" or "width" not {ip_axis}')

    if (ip_axis == 'length' and len_dominant) or (ip_axis == 'width' and not len_dominant):
        x_axis = True  # Direction of l
    else:
        x_axis = False  # Direction of b
    v = sl.poissons_ratio
    l_o_b = l / b
    if fd.depth:
        d_o_b = fd.depth / b
        n_emb = (1 + (0.33 + 1.34 / (1 + l_o_b)) * d_o_b ** 0.8)  # Note: emb_x = emb_y
    else:
        n_emb = 1.0
    if x_axis:  # x_axis
        k_h_surf = sl.g_mod * b / (2 - v) * (6.8 * l_o_b ** 0.65 + 2.4)
    else:
        k_h_surf = sl.g_mod * b / (2 - v) * (6.8 * l_o_b ** 0.65 + 0.8 * l_o_b + 1.6)
    f_dyn = 1.0  # Table 2-3a
    return k_h_surf * n_emb * f_dyn

