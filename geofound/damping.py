from geofound import tables_of_dyn_coefficients as tdc
import numpy as np

_pi = 3.14159265359


def calc_vert_via_gazetas_1991(sl, fd, a0, f_contact=1.0, saturated=False):
    v_s = sl.get_shear_vel(saturated=saturated)
    v_la = 3.4 / (_pi * (1 - sl.poissons_ratio)) * v_s
    if saturated:
        rho = sl.unit_sat_mass
    else:
        rho = sl.unit_dry_mass
    l = max([fd.length, fd.width]) / 2
    b = max([fd.length, fd.width]) / 2

    f_dyn = tdc.get_cz_gazetas_v_lte_0p4(a0, l / b)
    if sl.poissons_ratio > 0.4:
        czf = tdc.get_czf_gazetas_v_gt_0p4(a0, l / b)
        f_dyn *= czf

    if fd.depth:
        h = min([fd.height, fd.depth])
        a_w = 2 * h * (fd.width + fd.length) * f_contact
        c_emb = rho * v_s * a_w
    else:
        c_emb = 0.0
    return rho * v_la * fd.area * f_dyn + c_emb


def calc_vert_strip_via_gazetas_1991(sl, fd, a0, ip_axis='width', f_contact=1.0, saturated=False):
    v_s = sl.get_shear_vel(saturated=saturated)
    v_la = 3.4 / (_pi * (1 - sl.poissons_ratio)) * v_s
    if saturated:
        rho = sl.unit_sat_mass
    else:
        rho = sl.unit_dry_mass

    l_oop = 1.0
    l_ip = getattr(fd, ip_axis)
    if a0:
        f_dyn = tdc.get_cz_gazetas_v_lte_0p4(a0, lob=1000)
        if sl.poissons_ratio > 0.4:
            czf = tdc.get_czf_gazetas_v_gt_0p4(a0, lob=1000)
            f_dyn *= czf
    else:
        f_dyn = 1.0
    if fd.depth:
        h = min([fd.height, fd.depth])
        a_w = 2 * h * (l_ip + l_oop) * f_contact
        c_emb = rho * v_s * a_w
    else:
        c_emb = 0.0
    return rho * v_la * (l_ip * l_oop) * f_dyn + c_emb


def calc_horz_via_gazetas_1991(sl, fd, a0, ip_axis='width', f_contact=1.0, saturated=False):
    if saturated:
        rho = sl.unit_sat_mass
    else:
        rho = sl.unit_dry_mass
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
    if y_axis:
        f_dyn_v3 = tdc.get_cy_gazetas_v_e_0p3(a0, l / b)
        f_dyn_v5 = tdc.get_cy_gazetas_v_e_0p3(a0, l / b)
        f_dyn = np.interp(sl.poissons_ratio, [0.3, 0.5], [f_dyn_v3, f_dyn_v5])
    else:
        f_dyn = 1.0  # no dynamic effect
    if fd.depth:
        v_s = sl.get_shear_vel(saturated=saturated)
        v_la = 3.4 / (_pi * (1 - sl.poissons_ratio)) * v_s
        l_ip = getattr(fd, ip_axis)
        if ip_axis == 'width':
            l_oop = fd.length
        else:
            l_oop = fd.width
        h = min([fd.height, fd.depth])
        a_wc = 2 * h * l_oop * f_contact
        a_ws = 2 * h * l_ip * f_contact
        c_emb = rho * v_s * a_ws + rho * v_la * a_wc
    else:
        c_emb = 0.0
    return rho * sl.get_shear_vel(saturated=saturated) * fd.area * f_dyn + c_emb


def calc_horz_strip_via_gazetas_1991(sl, fd, a0, ip_axis='width', f_contact=1.0, saturated=False):
    if saturated:
        rho = sl.unit_sat_mass
    else:
        rho = sl.unit_dry_mass

    f_dyn_v3 = tdc.get_cy_gazetas_v_e_0p3(a0, lob=1000)
    f_dyn_v5 = tdc.get_cy_gazetas_v_e_0p3(a0, lob=1000)
    f_dyn = np.interp(sl.poissons_ratio, [0.3, 0.5], [f_dyn_v3, f_dyn_v5])
    l_oop = 1.0
    l_ip = getattr(fd, ip_axis)
    if fd.depth:
        v_s = sl.get_shear_vel(saturated=saturated)
        v_la = 3.4 / (_pi * (1 - sl.poissons_ratio)) * v_s
        h = min([fd.height, fd.depth])
        a_wc = 2 * h * l_oop * f_contact
        a_ws = 2 * h * l_ip * f_contact
        c_emb = rho * v_s * a_ws + rho * v_la * a_wc
    else:
        c_emb = 0.0
    return rho * sl.get_shear_vel(saturated=saturated) * (l_ip * l_oop) * f_dyn + c_emb


def calc_rot_via_gazetas_1991(sl, fd, a0, ip_axis='width', saturated=False, f_contact=1.0):
    v_la = 3.4 / (_pi * (1 - sl.poissons_ratio)) * sl.get_shear_vel(saturated=saturated)
    if saturated:
        rho = sl.unit_sat_mass
    else:
        rho = sl.unit_dry_mass
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
    if xx_axis:
        c_static = rho * v_la * i_bx
        f_dyn = tdc.get_crx_gazetas(a0, l / b)
        c_emb = 0.0
    else:
        c_static = rho * v_la * i_by
        f_dyn = tdc.get_cry_gazetas(a0, l / b)
        c_emb = 0.0  # TODO: this is wrong - but formula is very hard to interpret
    return c_static * f_dyn + c_emb


def calc_rot_strip_via_gazetas_1991(sl, fd, a0, ip_axis='width', saturated=False, f_contact=1.0):
    v_la = 3.4 / (_pi * (1 - sl.poissons_ratio)) * sl.get_shear_vel(saturated=saturated)
    if saturated:
        rho = sl.unit_sat_mass
    else:
        rho = sl.unit_dry_mass
    l_oop = 1.0
    l_ip = getattr(fd, ip_axis)
    i_bx = l_oop * l_ip ** 3 / 12

    c_static = rho * v_la * i_bx
    f_dyn = tdc.get_crx_gazetas(a0, lob=1000)
    c_emb = 0.0

    return c_static * f_dyn + c_emb


def calc_tors_via_gazetas_1991(sl, fd, a0, saturated=False):
    j_t = fd.i_ll + fd.i_ww
    if saturated:
        rho = sl.unit_sat_mass
    else:
        rho = sl.unit_dry_mass
    l = max([fd.length, fd.width]) / 2
    b = max([fd.length, fd.width]) / 2
    if a0:
        f_dyn = tdc.get_ct_gazetas(a0, l / b)
    else:
        f_dyn = 1.0
    return rho * sl.get_shear_vel(saturated=saturated) * j_t * f_dyn


def show_example():
    import geofound as gf

    sl = gf.create_soil()
    sl.g_mod = 30e6
    sl.poissons_ratio = 0.3
    sl.unit_dry_weight = 16.0e3

    lens = [0.5, 1.5, 10]
    length = 1.5
    width = 1.5
    depth = 0.0
    fd = gf.create_foundation(length, width, depth)
    a0 = 1
    c_rot = calc_rot_via_gazetas_1991(sl, fd, ip_axis='width', a0=a0)
    c_strip = calc_rot_strip_via_gazetas_1991(sl, fd, ip_axis='width', a0=a0)
    print(c_strip / c_rot * fd.length)
    c_vert = calc_vert_via_gazetas_1991(sl, fd, a0=a0)
    c_strip = calc_vert_strip_via_gazetas_1991(sl, fd, ip_axis='width', a0=a0)
    print(c_strip / c_vert * fd.length)
    c_horz = calc_horz_via_gazetas_1991(sl, fd, ip_axis='width', a0=a0)
    c_strip = calc_horz_strip_via_gazetas_1991(sl, fd, ip_axis='width', a0=a0)
    print(c_strip / c_horz * fd.length)
    calc_rot_via_gazetas_1991(sl, fd, ip_axis='width', a0=0.5)


if __name__ == '__main__':
    show_example()