import geofound as gf


def rotational_stiffness(sl, fd, ip_axis="width", axis=None, a0=0.0, method='gazetas_1991', **kwargs):
    """
    Rotation stiffness of foundation.

    Parameters
    ----------
    fd: Foundation object
    sl: Soil Object.
    ip_axis: str
        The axis that is in the plane of deformation
    axis: str
        The axis which it should be computed around (if not None, then ip_axis is ignored)

    k_f: float
        Rotational stiffness of the foundation
    """
    if method == 'gazetas_1991':
        return calc_rotational_via_gazetas_1991(sl, fd, ip_axis=ip_axis, axis=axis, a0=a0, **kwargs)
    else:
        return calc_rotational_via_pais_1988(sl, fd, ip_axis=ip_axis, axis=axis, a0=a0, **kwargs)


def calc_rotational_via_pais_1988(sl, fd, ip_axis='width', axis=None, a0=0.0, **kwargs):
    """
    Rotation stiffness of foundation from Pais and Kausel (1988).

    Parameters
    ----------
    fd: Foundation object
    sl: Soil Object
    ip_axis: str
        The axis that is in the plane of deformation
    axis: str
     The axis which it should be computed around (if not None, then ip_axis is ignored)
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
    if axis is not None:
        if axis == 'length':
            ip_axis = 'width'
        else:
            ip_axis = 'length'

    if fd.i_ll >= fd.i_ww:
        len_dominant = True
        l = fd.length * 0.5
        b = fd.width * 0.5
    else:
        len_dominant = False
        l = fd.width * 0.5
        b = fd.length * 0.5
    if (ip_axis == 'width' and len_dominant) or (ip_axis == 'length' and not len_dominant):
        x_axis = True
    else:
        x_axis = False

    v = sl.poissons_ratio
    n_emb = 1
    if x_axis:
        k_rx = 1 - 0.2 * a0
        # x-axis
        k_f_0 = (sl.g_mod * b ** 3 / (1 - v) * (3.2 * (l / b) + 0.8)) * k_rx
        if fd.depth > 0.0:
            n_emb = 1.0 + fd.depth / b + (1.6 / (0.35 + (l / b)) * (fd.depth / b) ** 2)
    else:
        k_ry = 1 - 0.3 * a0
        k_f_0 = (sl.g_mod * b ** 3 / (1 - v) * (3.73 * (l / b) ** 2.4 + 0.27)) * k_ry
        if fd.depth > 0.0:
            n_emb = 1.0 + fd.depth / b + (1.6 / (0.35 + (l / b) ** 4) * (fd.depth / b) ** 2)

    return k_f_0 * n_emb


def calc_rotational_via_gazetas_1991(sl, fd, ip_axis='width', axis="length", a0=0.0, f_contact=1.0, **kwargs):
    """
    Rotation stiffness of foundation from Gazetas (1991) and Mylonakis et al. (2006)

    Parameters
    ----------
    fd: Foundation object
    sl: Soil Object
    ip_axis: str
        The axis that is in the plane of deformation
    axis: str
     The axis which it should be computed around (if not None, then ip_axis is ignored)
    a0: float
        dynamic factor
    f_contact: float
        Effective sidewall contact scaling factor

    Returns
    -------
    k_f: float
        Rotational stiffness of the foundation
    """
    if not kwargs.get("disable_requires", False):
        gf.models.check_required(sl, ["g_mod", "poissons_ratio"])
        gf.models.check_required(fd, ["length", "width", "depth"])
    if axis is not None:
        if axis == 'length':
            ip_axis = 'width'
        else:
            ip_axis = 'length'

    if fd.i_ww >= fd.i_ll:
        len_dominant = True
        l = fd.length * 0.5
        b = fd.width * 0.5
        i_bx = fd.i_ll
        i_by = fd.i_ww
    else:
        len_dominant = False
        l = fd.width * 0.5
        b = fd.length * 0.5
        i_by = fd.i_ll
        i_bx = fd.i_ww
    if (ip_axis == 'width' and len_dominant) or (ip_axis == 'length' and not len_dominant):
        x_axis = True
    else:
        x_axis = False

    v = sl.poissons_ratio
    n_emb = 1.0
    if x_axis:
        k_rx = 1 - 0.2 * a0
        k_f_0 = (sl.g_mod / (1 - v) * i_bx ** 0.75 * (l / b) ** 0.25 * (2.4 + 0.5 * (b / l))) * k_rx
        if fd.depth > 0.0:
            dw = min(fd.height, fd.depth) * f_contact
            n_emb = 1 + 1.26 * (dw / b) ** 0.6 * (1.5 + (dw / fd.depth) ** 1.9 * (b / l) ** -0.6)
    else:
        k_ry = 1 - 0.3 * a0
        k_f_0 = (sl.g_mod / (1 - v) * i_by ** 0.75 * (3 * (l / b) ** 0.15)) * k_ry
        if fd.depth > 0.0:
            dw = min(fd.height, fd.depth) * f_contact
            n_emb = 1 + 0.92 * (dw / b) ** 0.6 * (1.0 + (dw / fd.depth) ** -0.2 * (b / l) ** 0.5)
    return k_f_0 * n_emb


def shear_stiffness(f_length, f_breadth, soil_g, soil_v):
    l = f_length * 0.5
    b = f_breadth * 0.5
    k_y = 2 * soil_g * l / (2 - soil_v) * (2.0 + 2.5 * (b / l) ** 0.85)
    k_shear = (k_y - (0.2 * soil_g * l) / (0.75 - soil_v) * (1.0 - b / l))
    return k_shear


def calc_shear_via_gazetas_1991(sl, fd, ip_axis='width', axis=None, a0=0.0, f_contact=1.0):
    """
    Calculate the shear stiffness for translation along an axis.

    Parameters
    ----------
    fd
    sl
    ip_axis: str
        The axis that is in the plane of deformation
    a0

    Returns
    -------

    """
    if axis is not None:
        ip_axis = axis
    # TODO: ADD depth correction
    n_emb = 1.0
    if fd.length >= fd.width:
        len_dominant = True
        l = fd.length * 0.5
        b = fd.width * 0.5
    else:
        len_dominant = False
        l = fd.width * 0.5
        b = fd.length * 0.5
    if (ip_axis == 'length' and len_dominant) or (ip_axis == 'width' and not len_dominant):
        y_axis = True  # Direction of l
    else:
        y_axis = False  # Direction of b
    v = sl.poissons_ratio
    k_y = 2 * sl.g_mod * l / (2 - v) * (2.0 + 2.5 * (b / l) ** 0.85)
    if y_axis is False:  # x_axis
        k_shear = (k_y - (0.2 * sl.g_mod * l) / (0.75 - v) * (1.0 - b / l))
        if fd.depth:
            z_w = (fd.depth + (fd.depth - fd.height)) / 2
            a_w = 2 * fd.height * (fd.width + fd.length) * f_contact
            n_emb = (1 + 0.15 * (fd.depth / l) ** 0.5) * (1 + 0.52 * (z_w * a_w / (l * b ** 2)) ** 0.4)
    else:
        if fd.depth:
            z_w = (fd.depth + (fd.depth - fd.height)) / 2
            a_w = 2 * fd.height * (fd.width + fd.length) * f_contact
            n_emb = (1 + 0.15 * (fd.depth / b) ** 0.5) * (1 + 0.52 * (z_w * a_w / (b * l ** 2)) ** 0.4)
        k_shear = k_y
    return k_shear * n_emb


def dep_rotational_stiffness(f_length, f_breadth, soil_g, soil_v):
    l = f_length * 0.5
    b = f_breadth * 0.5
    Iy = f_length ** 3 * f_breadth / 12
    k_f_0 = (soil_g / (1.0 - soil_v) * Iy ** 0.75 * (3 * (l / b) ** 0.15))
    return k_f_0


def get_vert_pais_1988(sl, fd, a0):
    v = sl.poissons_ratio
    l = fd.length * 0.5
    b = fd.width * 0.5
    k_v_0 = sl.g_mod * b / (1 - v) * (3.1 * (l / b) ** 0.75 + 1.6)
    return k_v_0


def calc_vert_via_gazetas_1991(sl, fd, a0=0.0, f_contact=1.0):
    v = sl.poissons_ratio
    l = fd.length * 0.5
    b = fd.width * 0.5
    k_v_0 = 2 * sl.g_mod * l / (1 - v) * (0.73 + 1.54 * (b / l) ** 0.75)
    n_emb = 1.0
    if fd.depth:
        a_w = 2 * fd.height * (fd.width + fd.length) * f_contact
        n_emb = (1 + fd.depth / (21 * b) * (1 + 1.3 * b / l)) * (1 + 0.2 * (a_w / (4 * b * l)) ** (2. / 3))
    return k_v_0 * n_emb


def calc_vert_via_pais_1988(sl, fd, a0=0.0):
    v = sl.poissons_ratio
    l = fd.length * 0.5
    b = fd.width * 0.5
    k_v_0 = sl.g_mod * b / (1 - v) * (3.1 * (l / b) ** 0.75 + 1.6)
    return k_v_0


def get_vert_gazetas_1991(sl, fd, a0):
    return calc_vert_via_gazetas_1991(sl, fd, a0)


def get_rot_via_gazetas_1991(sl, fd, ip_axis="length", a0=0.0, f_contact=1.0, **kwargs):
    return calc_rotational_via_gazetas_1991(sl, fd, ip_axis=ip_axis, a0=a0, f_contact=f_contact, **kwargs)
