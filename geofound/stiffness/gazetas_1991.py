import geofound as gf
from geofound import tables_of_dyn_coefficients as tdc


_pi = 3.14159265359


def calc_rot_strip_via_gazetas_1991(sl, fd, ip_axis=None, a0=0.0, f_contact=1.0, h_rigid=None, **kwargs):
    """
    Assumes out-of-plane is infinite, stiffness is returned as a per metre length

    Parameters
    ----------
    sl: Soil Object
    fd: Foundation object
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

    """
    if ip_axis is None:
        ip_axis = fd.ip_axis
    l_ip = getattr(fd, ip_axis)
    b = l_ip / 2
    # Note that Gazetas 1991 is x2 this, however, inconsistent with Gazetas (1983) table 4 pg 22
    # This formula is also consistent with num modelling
    k_strip = _pi * sl.g_mod * b ** 2 / (2 * (1 - sl.poissons_ratio))
    f_dyn = 1 - 0.2 * a0
    if fd.depth is not None and fd.depth != 0.0:
        if fd.depth < 0.0:
            raise ValueError(f'foundation depth must be zero or greater, not {fd.depth}')
        if f_contact == 0.0:
            dw = 0.0
        else:
            dw = min(fd.height, fd.depth) * f_contact
        # note: in ATC-40 the 1.26 factor is 2.52
        n_emb = 1 + 1.26 * (dw / b)
    else:
        n_emb = 1.
    if h_rigid:
        n_rigid = 1 + 0.2 * b / h_rigid
    else:
        n_rigid = 1
    return k_strip * n_emb * f_dyn * n_rigid


def calc_horz_via_gazetas_1991(sl, fd, ip_axis=None, a0=0.0, f_contact=1.0):
    """
    Calculate the shear stiffness for translation along an axis.

    Parameters
    ----------
    sl: Soil Object
    fd: Foundation object
    ip_axis: str
        The axis that is in the plane of deformation
    a0: float
        dynamic factor
    f_contact: float
        Effective sidewall contact scaling factor

    Returns
    -------

    """
    if ip_axis is None:
        ip_axis = fd.ip_axis
    n_emb = 1.0
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
    k_y = 2 * sl.g_mod * l / (2 - v) * (2.0 + 2.5 * (b / l) ** 0.85)
    if x_axis:  # x_axis
        k_shear = (k_y - (0.2 * sl.g_mod * l) / (0.75 - v) * (1.0 - b / l))
        f_dyn = 1.0
        if fd.depth:
            d_con = min([fd.height, fd.depth]) * f_contact
            z_w = fd.depth - d_con / 2
            a_w = 2 * d_con * (fd.width + fd.length)
            n_emb = (1 + 0.15 * (fd.depth / b) ** 0.5) * (1 + 0.52 * (z_w * a_w / (l * b ** 2)) ** 0.4)
            # TODO: add f_dyn for embedded case
    else:
        if a0:
            f_dyn = tdc.get_ky_gazetas(a0, l / b)
        else:
            f_dyn = 1.0
        if fd.depth:
            d_con = min([fd.height, fd.depth]) * f_contact
            z_w = fd.depth - d_con / 2
            a_w = 2 * d_con * (fd.width + fd.length)
            n_emb = (1 + 0.15 * (fd.depth / b) ** 0.5) * (1 + 0.52 * (z_w * a_w / (b * l ** 2)) ** 0.4)
            # TODO: add f_dyn for embedded case
        k_shear = k_y
    return k_shear * n_emb * f_dyn


def calc_shear_via_gazetas_1991(sl, fd, ip_axis=None, a0=0.0, f_contact=1.0):
    return calc_horz_via_gazetas_1991(sl, fd, ip_axis=ip_axis, a0=a0, f_contact=f_contact)


def calc_horz_strip_via_gazetas_1991(sl, fd, ip_axis=None, a0=0.0, f_contact=1.0, h_rigid=None):
    """
    Assumes out-of-plane is infinite, stiffness is returned as a per metre length

    Parameters
    ----------
    sl: Soil Object
    fd: Foundation object
    ip_axis: str
        The axis that is in the plane of deformation
    a0: float
        dynamic factor
    f_contact: float
        Effective sidewall contact scaling factor

    Returns
    -------

    """
    if ip_axis is None:
        ip_axis = fd.ip_axis
    l_ip = getattr(fd, ip_axis)
    b = l_ip / 2
    # note that Gazetas (1983) uses a factor of 2.1 instead of 2 - see pg 22
    k_h_0_strip = 2.0 * sl.g_mod / (2 - sl.poissons_ratio)

    if a0:
        f_dyn = tdc.get_ky_gazetas(a0, 1000)
    else:
        f_dyn = 1.0
    if fd.depth is not None and fd.depth != 0.0:
        if fd.depth < 0.0:
            raise ValueError(f'foundation depth must be zero or greater, not {fd.depth}')
        if f_contact == 0.0:
            dw = 0.0
        else:
            dw = min(fd.height, fd.depth) * f_contact
        n_emb = 1 + 1.26 * (dw / b)
    else:
        n_emb = 1.
    if h_rigid:
        n_rigid = 1 + 2.0 * b / h_rigid
    else:
        n_rigid = 1
    return k_h_0_strip * n_emb * f_dyn * n_rigid


def calc_vert_via_gazetas_1991(sl, fd, a0=None, f_contact=1.0, h_rigid=None):
    """
    Vertical stiffness of foundation

    Parameters
    ----------
    sl: Soil Object
    fd: Foundation object
    a0: float
        Dynamic factor
    f_contact: float
        Effective sidewall contact scaling factor

    Returns
    -------

    """
    v = sl.poissons_ratio
    l = max([fd.length * 0.5, fd.width * 0.5])
    b = min([fd.length * 0.5, fd.width * 0.5])
    k_n_0 = 2 * sl.g_mod * l / (1 - v) * (0.73 + 1.54 * (b / l) ** 0.75)

    if a0:
        if sl.poissons_ratio <= 0.4:
            f_dyn_surf = tdc.get_kz_gazetas_v_lte_0p4(a0, l / b)
        else:
            f_dyn_surf = tdc.get_kz_gazetas_v_gt_0p4(a0, l / b)
    else:
        f_dyn_surf = 1.
    if fd.depth is not None and fd.depth != 0.0:
        if fd.depth < 0.0:
            raise ValueError(f'foundation depth must be zero or greater, not {fd.depth}')
        if f_contact == 0.0:
            dw = 0.0
        else:
            dw = min([fd.height, fd.depth]) * f_contact
        a_w = 2 * dw * (fd.width + fd.length)
        chi = b / l
        n_emb = (1 + fd.depth / (21 * b) * (1 + 1.3 * chi)) * (1 + 0.2 * (a_w / (4 * b * l)) ** (2. / 3))
        if a0:
            h = min([fd.height, fd.depth])
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
            f_dyn_emb = 1
    else:
        n_emb = 1.0
        f_dyn_emb = 1.0
    if h_rigid:
        n_rigid = 1 + (b / h_rigid) / (0.5 + b / l)
    else:
        n_rigid = 1

    return k_n_0 * n_emb * f_dyn_surf * f_dyn_emb * n_rigid


def calc_vert_strip_via_gazetas_1991(sl, fd, ip_axis=None, a0=0.0, f_contact=1.0, h_rigid=None, **kwargs):
    """
    Vertical stiffness per metre of a footing with infinite out-of-plane length

    Parameters
    ----------
    sl: Soil Object
    fd: Foundation object
    a0: float
        Dynamic factor
    f_contact: float
        Effective sidewall contact scaling factor
    f_contact
    h_rigid
    kwargs

    Returns
    -------

    """
    if ip_axis is None:
        ip_axis = fd.ip_axis
    l_ip = getattr(fd, ip_axis)
    b = l_ip / 2
    v = sl.poissons_ratio
    # k_n_0_strip = 2 * 0.73 * sl.g_mod / (1 - v)  # from Gazetas (1991) table 15.3
    k_n_0_strip = 1.23 * sl.g_mod / (1 - v)  # from Gazetas (1983) table 4 (pg 22)
    # note that the factor 1.23 appears to be consistent with FEM modelling
    if a0:
        lob = 1000
        if sl.poissons_ratio <= 0.4:
            f_dyn_surf = tdc.get_kz_gazetas_v_lte_0p4(a0, lob)
        else:
            f_dyn_surf = tdc.get_kz_gazetas_v_gt_0p4(a0, lob)
    else:
        f_dyn_surf = 1.
    if fd.depth is not None and fd.depth != 0.0:
        if fd.depth < 0.0:
            raise ValueError(f'foundation depth must be zero or greater, not {fd.depth}')
        if f_contact == 0.0:
            dw = 0.0
        else:
            dw = min([fd.height, fd.depth]) * f_contact
        a_w = 2 * dw * fd.width
        n_emb = (1 + fd.depth / (21 * b)) * (1 + 0.2 * (a_w / (2 * b)) ** (2. / 3))
        if v <= 0.4:
            f_dyn_full_emb = 1 - 0.09 * (fd.depth / b) ** (3. / 4) * a0 ** 2
            f_dyn_trench = 1 + 0.09 * (fd.depth / b) ** (3. / 4) * a0 ** 2
        else:
            f_dyn_full_emb = 1 - 0.35 * (fd.depth / b) ** 0.5 * a0 ** 3.5
            f_dyn_trench = f_dyn_surf  # Assume the same as the surface - this is currently not provided
            f_dyn_surf = 1.
        # interpolate between
        f_dyn_emb = f_dyn_trench + (f_dyn_full_emb - f_dyn_trench) * (dw / fd.depth)
    else:
        n_emb = 1
        f_dyn_emb = 1.0
    if h_rigid:
        n_rigid = 1 + 3.5 * b / h_rigid
    else:
        n_rigid = 1
    return k_n_0_strip * n_emb * f_dyn_surf * f_dyn_emb * n_rigid


def get_vert_gazetas_1991(sl, fd, a0):
    gf.exceptions.deprecation('get_vert_gazetas_1991')
    return calc_vert_via_gazetas_1991(sl, fd, a0)


def get_rot_via_gazetas_1991(sl, fd, ip_axis=None, a0=0.0, f_contact=1.0, **kwargs):
    gf.exceptions.deprecation('get_rot_via_gazetas_1991')
    return calc_rot_via_gazetas_1991(sl, fd, ip_axis=ip_axis, a0=a0, f_contact=f_contact, **kwargs)

def calc_rotational_via_gazetas_1991(sl, fd, ip_axis=None, a0=0.0, f_contact=1.0, **kwargs):
    return calc_rot_via_gazetas_1991(sl, fd, ip_axis=ip_axis, a0=a0, f_contact=f_contact, **kwargs)


def calc_norm_stiff_of_rect_via_gazetas_1991(loop_o_lip):
    """
    Calculate the normalised rotational to vertical stiffness ratio

    norm_k = K_M / (K_N * L_ip ** 2)

    Parameters
    ----------
    loop_o_lip: float
        Length in-plane divided by length out-of-plane
    """
    length = 1
    width = loop_o_lip
    depth = 0.0
    sl = gf.create_soil()
    sl.g_mod = 30e6
    sl.poissons_ratio = 0.3
    fd = gf.create_foundation(length, width, depth)
    ip_axis = 'length'
    k_rot = calc_rot_via_gazetas_1991(sl, fd, ip_axis)
    k_vert = calc_vert_via_gazetas_1991(sl, fd)
    n_rat = k_rot / k_vert / fd.length ** 2
    return n_rat
    # return np.where(loop_o_lip < 1,
    #                 2.48 / loop_o_lip * (1 + 5 * loop_o_lip) / (0.73 + 1.55 / loop_o_lip ** 0.75),
    #                 14.88 * loop_o_lip ** 0.6 / (0.73 + 1.55 * loop_o_lip ** 0.75))


def calc_rot_via_gazetas_1991(sl, fd, ip_axis=None, a0=0.0, f_contact=1.0, **kwargs):
    """
    Rotation stiffness of foundation from Gazetas (1991) and Mylonakis et al. (2006)

    Parameters
    ----------
    sl: Soil Object
    fd: Foundation object
    ip_axis: str
        The axis that is in the plane of deformation
    a0: float
        dynamic factor
    f_contact: float
        Effective sidewall contact scaling factor

    Returns
    -------
    k_f: float
        Rotational stiffness of the foundation
    """
    h_rigid = kwargs.get("h_rigid", None)
    if not kwargs.get("disable_requires", False):
        gf.models.check_required(sl, ["g_mod", "poissons_ratio"])
        gf.models.check_required(fd, ["length", "width", "depth"])
    if ip_axis is None:
        ip_axis = fd.ip_axis
    if ip_axis is None:
        raise ValueError('ip_axis must be set')

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
    if ip_axis is None:
        ip_axis = fd.ip_axis
    if ip_axis not in ['length', 'width']:
        raise ValueError(f'Must set ip_axis to either "length" or "width" not {ip_axis}')
    if (ip_axis == 'width' and len_dominant) or (ip_axis == 'length' and not len_dominant):
        xx_axis = True  # weaker rotation (rotation about longest length)
    else:
        xx_axis = False

    v = sl.poissons_ratio
    n_emb = 1.0
    if xx_axis:
        f_dyn = 1 - 0.2 * a0
        k_static_surf = (sl.g_mod / (1 - v) * i_bx ** 0.75 * (l / b) ** 0.25 * (2.4 + 0.5 * (b / l)))
        if fd.depth is not None and fd.depth != 0.0:
            if fd.depth < 0.0:
                raise ValueError(f'foundation depth must be zero or greater, not {fd.depth}')
            if f_contact == 0.0:
                dw = 0.0
            else:
                dw = min(fd.height, fd.depth) * f_contact
            if dw == 0:
                n_emb = 1
            else:
                # note: in ATC-40 the 1.26 factor is 2.52
                n_emb = 1 + 1.26 * (dw / b) * (1. + (dw / b) * (dw / fd.depth) ** -0.2 * (b / l) ** 0.5)
    else:  # yy_axis (rotation about y-axis)
        if v < 0.45:
            f_dyn = 1 - 0.3 * a0  # Note the f_dyn_emb = 1.0
        else:
            f_dyn = 1 - 0.25 * a0 * (l / b) ** 0.3
        k_static_surf = (sl.g_mod / (1 - v) * i_by ** 0.75 * (3 * (l / b) ** 0.15))
        if fd.depth is not None and fd.depth != 0.0:
            if fd.depth < 0.0:
                raise ValueError(f'foundation depth must be zero or greater, not {fd.depth}')
            if f_contact == 0.0:
                dw = 0.0
            else:
                dw = min(fd.height, fd.depth) * f_contact
            if dw == 0.0:
                n_emb = 1
            else:
                # Note that the original Gazetas (1991) paper has (dw / L) ** 1.9 * (dw / D) ** -0.6, also in ATC-40
                # From Gazetas (1983) it is explained that strip has much less embedment effect that circular, since less sidewall
                n_emb = 1 + 0.92 * (dw / b) ** 0.6 * (1.5 + (dw / l) ** 1.9 * (dw / fd.depth) ** -0.6)
                # if fd.depth / b > 2. / 3:
                #     warnings.warn('D/B should be less than or equal to 2/3 - See Gazetas (1983) Table 9', EquationWarning, stacklevel=2)
                # whereas Mylonakis has this form:
                # n_emb = 1 + 0.92 * (dw / b) ** 0.6 * (1.5 + (dw / fd.depth) ** 1.9 * (b / l) ** -0.6)
    if h_rigid:
        n_rigid = 1 + 0.2 * b / h_rigid
    else:
        n_rigid = 1
    return k_static_surf * f_dyn * n_emb * n_rigid
