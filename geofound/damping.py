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
