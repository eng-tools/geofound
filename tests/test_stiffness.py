import geofound
import numpy as np


def test_rotation():
    length = 3
    width = 2
    depth = 0.0
    sl = geofound.create_soil()
    sl.g_mod = 30e6
    sl.poissons_ratio = 0.3
    fd = geofound.create_foundation(length, width, depth)
    i_ll = width ** 3 * length / 12
    i_ww = length ** 3 * width / 12
    assert geofound.isclose(fd.i_ll, i_ll, rel_tol=0.001)
    assert geofound.isclose(fd.i_ww, i_ww, rel_tol=0.001)
    assert geofound.isclose(geofound.stiffness.rotational_stiffness(sl, fd, axis="length"), 218027424.1324, rel_tol=0.01)
    assert geofound.isclose(geofound.stiffness.rotational_stiffness(sl, fd, axis="width"), 422150729.0333, rel_tol=0.01)


def test_stiffness_ratio_formula():
    length = 3
    width = 2
    depth = 0.0
    sl = geofound.create_soil()
    sl.g_mod = 30e6
    sl.poissons_ratio = 0.3
    fd = geofound.create_foundation(length, width, depth)
    ip_axis = 'length'
    k_rot = geofound.stiffness.calc_rot_via_gazetas_1991(sl, fd, ip_axis)
    k_vert = geofound.stiffness.calc_vert_via_gazetas_1991(sl, fd)
    n_rat = k_rot / k_vert / fd.length ** 2
    loop_o_lip = fd.width / fd.length
    loop_o_lip = 1. / loop_o_lip
    k1 = 1.24 / loop_o_lip * (1 + 5 * loop_o_lip) / (0.73 + 1.55 / loop_o_lip ** 0.75)
    k2 = 7.44 * loop_o_lip ** 0.6 / (0.73 + 1.55 * loop_o_lip ** 0.75)
    n_rat1 = geofound.stiffness.calc_norm_stiff_of_rect_via_gazetas_1991(loop_o_lip=fd.width / fd.length)
    assert np.isclose(n_rat, n_rat1), (n_rat, n_rat1)


if __name__ == '__main__':
    test_stiffness_ratio_formula()
