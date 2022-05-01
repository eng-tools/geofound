import numpy as np

import geofound as gf


def calc_rot_via_pais_1988(sl, fd, ip_axis='width', a0=0.0, **kwargs):
    """
    Rotation stiffness of foundation from Pais and Kausel (1988).

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
    n_emb = 1
    if xx_axis:
        # x-axis
        k_rx_dyn_surf = 1 - ((0.55 + 0.01 * (l / b - 1) ** 0.5) * a0 ** 2 / ((2.4 - 0.4 / (l / b) ** 3) + a0 ** 2))
        k_f_0_static_surf = (sl.g_mod * b ** 3 / (1 - v) * (3.2 * (l / b) + 0.8))
        if fd.depth is not None and fd.depth != 0.0:
            d = fd.depth
            if fd.depth < 0.0:
                raise ValueError(f'foundation depth must be zero or greater, not {fd.depth}')
            n_emb = 1.0 + fd.depth / b + (1.6 / (0.35 + (l / b)) * (fd.depth / b) ** 2)
            chi = np.sqrt(2 * (1 - v) / (1 - 2 * v))
            k_f_0_static = k_f_0_static_surf * n_emb
            k_rx_dyn = ((4. / 3) * ((d / b) + (d/b)**3 + chi * (l / b) * (d / b) ** 3 + 3 * (d / b) * (l / b) + chi * (l / b)) * a0 ** 2 / (k_f_0_static / (sl.g_mod * b **3)) * ((1.8 / (1 + 1.75 * (l/b - 1)) + a0 **2)) +
                        (4. / 3) * (chi * l / b + 1) * (d/b)**3 / (k_f_0_static / (sl.g_mod * b **3))) * (a0 / (2 * k_rx_dyn_surf))
        else:
            k_f_0_static = k_f_0_static_surf
            k_rx_dyn = k_rx_dyn_surf
        k_f_0 = k_f_0_static * k_rx_dyn
    else:
        k_ry_dyn_surf = 1 - (0.55 * a0 ** 2 / ((0.6 + 1.4 / (l / b) ** 3) + a0 ** 2))
        k_f_0_static_surf = (sl.g_mod * b ** 3 / (1 - v) * (3.73 * (l / b) ** 2.4 + 0.27))
        if fd.depth is not None and fd.depth != 0.0:
            d = fd.depth
            if fd.depth < 0.0:
                raise ValueError(f'foundation depth must be zero or greater, not {fd.depth}')
            n_emb = 1.0 + d / b + (1.6 / (0.35 + (l / b) ** 4) * (d / b) ** 2)
            k_f_0_static = k_f_0_static_surf * n_emb
            k_ry_dyn = ((4. / 3) * ((l / b) ** 3 * (d / b) + chi * (l / b) * (d / b) ** 3 + (d/b)**3 + 3 * (d / b) * (l / b) ** 2 + chi * (
                            l / b) ** 3) * a0 ** 2 / (k_f_0_static / (sl.g_mod * b ** 3)) * ((1.8 / (1 + 1.75 * (l / b - 1)) + a0 ** 2)) +
                        (4. / 3) * (l / b + chi) * (d / b) ** 3 / (k_f_0_static / (sl.g_mod * b ** 3))) * (
                                   a0 / (2 * k_ry_dyn_surf))
        else:
            k_f_0_static = k_f_0_static_surf
            k_ry_dyn = k_ry_dyn_surf
        k_f_0 = k_f_0_static * k_ry_dyn

    return k_f_0 * n_emb


def calc_vert_via_pais_1988(sl, fd, a0=0):
    """
    Vertical stiffness of foundation from Pais and Kausel (1988).

    Parameters
    ----------
    sl
    fd
    a0

    Returns
    -------

    """
    v = sl.poissons_ratio
    l = fd.length * 0.5
    b = fd.width * 0.5
    k_v_0 = sl.g_mod * b / (1 - v) * (3.1 * (l / b) ** 0.75 + 1.6)
    kz = 1 - ((0.4 + 0.2 / (l / b)) * a0 ** 2 / (10 / (1 + 3 * (l / b) - 1) + a0 ** 2))
    n_emb = 1
    if fd.depth is not None and fd.depth != 0.0:
        if fd.depth < 0.0:
            raise ValueError(f'foundation depth must be zero or greater, not {fd.depth}')
        n_emb = 1 + (0.25 + 0.25 / (l / b)) * (fd.depth / b) ** 0.8

    if a0:
        f_dyn_surf = 1. - ((0.4 + 0.2 / (l / b) * a0 ** 2) / (10 / (1 + 3 * (l / b) - 1)) + a0 ** 2)

    else:
        f_dyn_surf = 1
    return k_v_0 * kz * n_emb * f_dyn_surf