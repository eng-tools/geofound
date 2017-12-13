
def rotational_stiffness(sp, fd, axis="length", a0=0.0):
    """
    Rotation stiffness of foundation.
    :param fd: Foundation object
    :param sp: Soil Object.
    :param axis: The axis which it should be computed around
    :return:
    """
    if fd.depth > 0.0:
        pass
    l = fd.length * 0.5
    b = fd.width * 0.5
    v = sp.poissons_ratio
    if axis == "length":
        i_bx = fd.i_ll
        k_rx = 1 - 0.2 * a0
        k_f_0 = (sp.g_mod / (1 - v) * i_bx ** 0.75 * (l / b) ** 0.25 * (2.4 + 0.5 * (b / l))) * k_rx
    else:
        i_by = fd.i_ww
        k_ry = 1 - 0.3 * a0
        k_f_0 = (sp.g_mod / (1 - v) * i_by ** 0.75 * (3 * (l / b) ** 0.15)) * k_ry
    return k_f_0


def shear_stiffness(f_length, f_breadth, soil_g, soil_v):
    l = f_length * 0.5
    b = f_breadth * 0.5
    k_y = 2 * soil_g * l / (2 - soil_v) * (2.0 + 2.5 * (b / l) ** 0.85)
    k_shear = (k_y - (0.2 * soil_g * l) / (0.75 - soil_v) * (1.0 - b / l))
    return k_shear


def dep_rotational_stiffness(f_length, f_breadth, soil_g, soil_v):
    l = f_length * 0.5
    b = f_breadth * 0.5
    Iy = f_length ** 3 * f_breadth / 12
    k_f_0 = (soil_g / (1.0 - soil_v) * Iy ** 0.75 * (3 * (l / b) ** 0.15))
    return k_f_0
