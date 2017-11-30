
def rotational_stiffness(sp, fd, axis="length"):
    """
    Rotation stiffness of mat foundation.
    :param fd: Foundation object
    :param sp: Soil Object.
    :param axis: The axis which it should be computed around
    :return:
    """
    if axis == "length":
        l = fd.length * 0.5
        b = fd.width * 0.5
        i_y = fd.i_ll
    else:
        l = fd.width * 0.5
        b = fd.length * 0.5
        i_y = fd.i_ww
    k_f_0 = (sp.g_mod / (1 - sp.poissons_ratio) * i_y ** 0.75 * (3 * (l / b) ** 0.15))
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
