import numpy as np

from geofound.output import log


def settlement_schmertmann(sp, fd, load, youngs_modulus_soil, **kwargs):
    """
    Calculates the settlement of a shallow foundation (Schmertmann, 19XX).

    :param sp: Soil Profile object
    :param fd: Foundation object
    :param load:
    :param youngs_modulus_soil: The Young's modulus of the soil.
    :param kwargs:
    :return: float, the settlement.
    """
    length = float(fd.length)
    breadth = float(fd.width)
    depth = float(fd.depth)
    load = float(load)

    sp.gwl = kwargs.get("gwl", sp.gwl)
    sp.unit_sat_weight = kwargs.get("unit_sat_weight", sp.unit_sat_weight)
    verbose = kwargs.get("verbose", 0)
    years = kwargs.get("years", 0)

    q = load / (length * breadth)
    sigma_v0_eff = (sp.unit_dry_weight * min(depth, sp.gwl) +
                    (sp.unit_sat_weight - 9.8) * max([0, depth - sp.gwl]))
    delta_q = q - sigma_v0_eff

    # EMBEDMENT FACTOR
    c_1 = max(1 - 0.5 * (sigma_v0_eff / delta_q), 0.5)
    # CREEP FACTOR
    if years == 0:
        c_2 = 1.0
    else:
        c_2 = 1.0 + 0.2 * np.log10(years / 0.1)
    # SHAPE FACTOR
    long = max(length, breadth)
    short = min(length, breadth)
    c_3 = max(1.03 - 0.03 * (long / short), 0.73)
    # Peak settlement index
    if long / short > 10:
        zp = short + depth
        z_top = 0.2
        z_bottom = 4 * short + depth
    else:
        z_top = 0.1
        zp = 0.5 * short + depth
        z_bottom = 2 * short + depth

    sigma_vp_eff = (sp.unit_dry_weight * min(zp, sp.gwl) +
                    (sp.unit_sat_weight - 9.8) * max([0, zp - sp.gwl]))
    i_zp = 0.5 + 0.1 * (delta_q / sigma_vp_eff) ** 0.5

    i_z_top = (i_zp + z_top) / 2
    i_z_bottom = i_zp / 2

    settlement = (c_1 * c_2 * c_3 * delta_q *
                  (i_z_top * (zp - depth) + i_z_bottom * (z_bottom - zp)) / youngs_modulus_soil)

    if verbose:
        log("delta_q:", delta_q)
        log("c_1:", c_1)
        log("c_2:", c_2)
        log("c_3:", c_3)
        log("zp:", zp)
        log("sigma_vp_eff:", sigma_vp_eff)
        log("i_zp:", i_zp)
        log("i_z_top:", i_z_top)
        log("i_z_bottom:", i_z_bottom)
        log("settlement:", settlement)
    return settlement
