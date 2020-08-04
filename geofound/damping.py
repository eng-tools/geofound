import numpy as np
from geofound import tables_of_dyn_coefficients as tdc


def calc_vert_via_gazetas_1991(sl, fd, a0=None, saturated=False):
    v_la = 3.4 / (np.pi * (1 - sl.poissons_ratio)) * sl.get_shear_vel(saturated=saturated)
    if saturated:
        rho = sl.unit_sat_mass
    else:
        rho = sl.unit_dry_mass
    l = max([fd.length, fd.width]) / 2
    b = max([fd.length, fd.width]) / 2
    if a0:
        f_dyn = tdc.get_cz_gazetas_v_lte_0p4(a0, l / b)
        if sl.poissons_ratio > 0.4:
            czf = tdc.get_czf_gazetas_v_gt_0p4(a0, l / b)
            f_dyn *= czf
    else:
        f_dyn = 1.0
    return rho * v_la * fd.area * f_dyn


def calc_horz_via_gazetas_1991(sl, fd, a0=None, ip_axis='width', saturated=False):
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
        if a0:
            f_dyn_v3 = tdc.get_cy_gazetas_v_e_0p3(a0, l / b)
            f_dyn_v5 = tdc.get_cy_gazetas_v_e_0p3(a0, l / b)
            f_dyn = np.interp(sl.poissons_ratio, [0.3, 0.5], f_dyn_v3, f_dyn_v5)
        else:
            f_dyn = 1.0
    else:
        f_dyn = 1.0  # no dynamic effective
    return rho * sl.get_shear_vel(saturated=saturated) * fd.area * f_dyn


def calc_rot_via_gazetas_1991(sl, fd, a0=None, ip_axis='width', saturated=False):
    v_la = 3.4 / (np.pi * (1 - sl.poissons_ratio)) * sl.get_shear_vel(saturated=saturated)
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
    else:
        c_static = rho * v_la * i_by
        f_dyn = tdc.get_cry_gazetas(a0, l / b)
    return c_static * f_dyn
