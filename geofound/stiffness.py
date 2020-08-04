import geofound as gf
from geofound.exceptions import EquationWarning
import warnings
from geofound import tables_of_dyn_coefficients as tdc


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

    if fd.i_ww >= fd.i_ll:
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
        k_rx = 1 - ((0.55 + 0.01 * (l / b - 1) ** 0.5) * a0 ** 2 / ((2.4 - 0.4 / (l / b) ** 3) + a0 ** 2))
        # x-axis
        k_f_0 = (sl.g_mod * b ** 3 / (1 - v) * (3.2 * (l / b) + 0.8)) * k_rx
        if fd.depth > 0.0:
            n_emb = 1.0 + fd.depth / b + (1.6 / (0.35 + (l / b)) * (fd.depth / b) ** 2)
    else:
        k_ry = 1 - (0.55 * a0 ** 2 / ((0.6 + 1.4 / (l / b) ** 3) + a0 ** 2))
        k_f_0 = (sl.g_mod * b ** 3 / (1 - v) * (3.73 * (l / b) ** 2.4 + 0.27)) * k_ry
        if fd.depth > 0.0:
            n_emb = 1.0 + fd.depth / b + (1.6 / (0.35 + (l / b) ** 4) * (fd.depth / b) ** 2)

    return k_f_0 * n_emb


def calc_rotational_via_gazetas_1991(sl, fd, ip_axis='width', axis=None, a0=0.0, f_contact=1.0, **kwargs):
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

    if fd.length >= fd.width:
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
        xx_axis = True  # weaker rotation
    else:
        xx_axis = False

    v = sl.poissons_ratio
    n_emb = 1.0
    if xx_axis:
        k_rx = 1 - 0.2 * a0
        k_f_0 = (sl.g_mod / (1 - v) * i_bx ** 0.75 * (l / b) ** 0.25 * (2.4 + 0.5 * (b / l))) * k_rx
        if fd.depth > 0.0:
            dw = min(fd.height, fd.depth) * f_contact
            #note: in ATC-40 the 1.26 factor is 2.52
            n_emb = 1 + 1.26 * (dw / b) * (1. + (dw / b) * (dw / fd.depth) ** -0.2 * (b / l) ** 0.5)
    else:  # yy_axis (rotation about y-axis)
        if v < 0.45:
            k_ry = 1 - 0.3 * a0
        else:
            k_ry = 1 - 0.25 * a0 * (l / b) ** 0.3
        k_f_0 = (sl.g_mod / (1 - v) * i_by ** 0.75 * (3 * (l / b) ** 0.15)) * k_ry
        if fd.depth > 0.0:
            dw = min(fd.height, fd.depth) * f_contact
            # Note that the original Gazetas (1991) paper has (dw / L) ** 1.9 * (dw / D) ** -0.6, also in ATC-40
            # From Gazetas (1983) it is explained that strip has much less embedment effect that circular, since less sidewall
            n_emb = 1 + 0.92 * (dw / b) ** 0.6 * (1.5 + (dw / l) ** 1.9 * (dw / fd.depth) ** -0.6)
            if fd.depth / b > 2. / 3:
                warnings.warn('D/B should be less than or equal to 2/3 - See Gazetas (1983)', EquationWarning, stacklevel=2)
            # whereas Mylonakis has this form:
            # n_emb = 1 + 0.92 * (dw / b) ** 0.6 * (1.5 + (dw / fd.depth) ** 1.9 * (b / l) ** -0.6)
    return k_f_0 * n_emb


def calc_rotational_strip_via_gazetas_1991(sl, fd, ip_axis='width', a0=0.0, f_contact=1.0, h_rigid=None, **kwargs):
    """
    Assumes out-of-plane is infinite, stiffness is returned as a per metre length

    Parameters
    ----------
    sl
    fd
    ip_axis
    axis
    a0
    f_contact
    kwargs

    Returns
    -------

    """
    l_ip = getattr(fd, ip_axis)
    b = l_ip / 2
    k_strip = 0.73 * sl.g_mod * b / (1 - sl.poissons_ratio)

    if 1:
        pass
    else:  # yy_axis (rotation about y-axis)
        if v < 0.45:
            k_ry = 1 - 0.3 * a0
        else:
            k_ry = 1 - 0.25 * a0 * (l / b) ** 0.3
        k_f_0 = (sl.g_mod / (1 - v) * i_by ** 0.75 * (3 * (l / b) ** 0.15)) * k_ry
        if fd.depth > 0.0:
            dw = min(fd.height, fd.depth) * f_contact
            # Note that the original Gazetas (1991) paper has (dw / L) ** 1.9 * (dw / D) ** -0.6, also in ATC-40
            # From Gazetas (1983) it is explained that strip has much less embedment effect that circular, since less sidewall
            n_emb = 1 + 0.92 * (dw / b) ** 0.6 * (1.5 + (dw / l) ** 1.9 * (dw / fd.depth) ** -0.6)
            if fd.depth / b > 2. / 3:
                warnings.warn('D/B should be less than or equal to 2/3 - See Gazetas (1983)', EquationWarning, stacklevel=2)
            # whereas Mylonakis has this form:
            # n_emb = 1 + 0.92 * (dw / b) ** 0.6 * (1.5 + (dw / fd.depth) ** 1.9 * (b / l) ** -0.6)
    return k_f_0 * n_emb


def shear_stiffness(f_length, f_breadth, soil_g, soil_v):
    gf.exceptions.deprecation('shear_stiffness')
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
    if a0:
        raise ValueError('dynamic stiffness not implemented')
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
            h = min([fd.height, fd.depth])
            a_w = 2 * h * (fd.width + fd.length) * f_contact
            n_emb = (1 + 0.15 * (fd.depth / l) ** 0.5) * (1 + 0.52 * (z_w * a_w / (l * b ** 2)) ** 0.4)
    else:
        if fd.depth:
            z_w = (fd.depth + (fd.depth - fd.height)) / 2
            h = min([fd.height, fd.depth])
            a_w = 2 * h * (fd.width + fd.length) * f_contact
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
    gf.exceptions.deprecation('get_vert_pais_1988 has deprecated. use calc_vert_via_pais_1988')
    v = sl.poissons_ratio
    l = fd.length * 0.5
    b = fd.width * 0.5
    k_v_0 = sl.g_mod * b / (1 - v) * (3.1 * (l / b) ** 0.75 + 1.6)
    return k_v_0


def calc_vert_via_pais_1988(sl, fd, a0=0):
    v = sl.poissons_ratio
    l = fd.length * 0.5
    b = fd.width * 0.5
    k_v_0 = sl.g_mod * b / (1 - v) * (3.1 * (l / b) ** 0.75 + 1.6)
    kz = 1 - ((0.4 + 0.2 / (l / b)) * a0 ** 2 / (10 / (1 + 3 * (l / b) - 1) + a0 ** 2))
    n_emb = 1
    if fd.depth:
        n_emb = 1 + (0.25 + 0.25 / (l / b)) * (fd.depth / b) ** 0.8
    return k_v_0 * kz * n_emb


def calc_vert_via_gazetas_1991(sl, fd, a0=None, f_contact=1.0):
    v = sl.poissons_ratio
    l = fd.length * 0.5
    b = fd.width * 0.5
    k_v_0 = 2 * sl.g_mod * l / (1 - v) * (0.73 + 1.54 * (b / l) ** 0.75)

    if a0:
        if sl.poissons_ratio <= 0.4:
            f_dyn_surf = tdc.get_kz_gazetas_v_lte_0p4(a0, l / b)
        else:
            f_dyn_surf = tdc.get_kz_gazetas_v_gt_0p4(a0, l / b)
    else:
        f_dyn_surf = 1.
    if fd.depth:
        h = min([fd.height, fd.depth])
        a_w = 2 * h * (fd.width + fd.length) * f_contact
        chi = b / l
        n_emb = (1 + fd.depth / (21 * b) * (1 + 1.3 * chi)) * (1 + 0.2 * (a_w / (4 * b * l)) ** (2. / 3))
        if v <= 0.4:
            f_dyn_full_emb = 1 - 0.09 * (fd.depth / b) ** (3. / 4) * a0 ** 2
            f_dyn_trench = 1 + 0.09 * (fd.depth / b) ** (3. / 4) * a0 ** 2
        else:

            if l / b < 2.5:
                f_dyn_full_emb = 1 - 0.09 * (fd.depth / b) ** (3. / 4) * a0 ** 2
            else:
                f_dyn_full_emb = 1 - 0.35 * (fd.depth / b) ** 0.5 * a0 ** 3.5
            f_dyn_trench = f_dyn_surf  # Assume the same as the surface - this is currently not provided
            f_dyn_surf = 1.
        # interpolate between
        f_dyn_emb = f_dyn_trench + (f_dyn_full_emb - f_dyn_trench) * ((h * f_contact) / fd.depth)
    else:
        n_emb = 1.0
        f_dyn_emb = 1.0

    return k_v_0 * n_emb * f_dyn_surf * f_dyn_emb


def calc_vert_strip_via_gazetas_1991(sl, fd, ip_axis='width', a0=0.0, f_contact=1.0, h_rigid=None, **kwargs):
    l_ip = getattr(fd, ip_axis)
    b = l_ip / 2
    k_strip = 0.73 * sl.g_mod * b / (1 - sl.poissons_ratio)
    if fd.depth:
        h = min([fd.height, fd.depth])
        a_w = 2 * fd.height * fd.width * f_contact
        n_emb = (1 + fd.depth / (21 * b)) * (1 + 0.2 * (a_w / (2 * b)) ** (2. / 3))
        f_dyn_full_emb = 1 - 0.09 * (fd.depth / b) ** (3. / 4) * a0 ** 2
        f_dyn_trench = 1 + 0.09 * (fd.depth / b) ** (3. / 4) * a0 ** 2
        # interpolate between
        f_dyn_full = f_dyn_trench + (f_dyn_full_emb - f_dyn_trench) * ((h * f_contact) / fd.depth)


    if fd.depth > 0.0:
        dw = min(fd.height, fd.depth) * f_contact
        if h_rigid is None:
            n_rigid = 1
        else:
            n_rigid = 1 + 3.5 * b / (h_rigid - fd.depth)
        n_emb = (1 + 0.2 * (dw / b) ** (2. / 3))


def get_vert_gazetas_1991(sl, fd, a0):
    gf.exceptions.deprecation('get_vert_gazetas_1991')
    return calc_vert_via_gazetas_1991(sl, fd, a0)


def get_rot_via_gazetas_1991(sl, fd, ip_axis="length", a0=0.0, f_contact=1.0, **kwargs):
    gf.exceptions.deprecation('get_rot_via_gazetas_1991')
    return calc_rotational_via_gazetas_1991(sl, fd, ip_axis=ip_axis, a0=a0, f_contact=f_contact, **kwargs)
