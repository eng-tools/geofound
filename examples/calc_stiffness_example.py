import geofound as gf


def show_example():
    length = 1.5
    width = 1.5
    depth = 0.0
    sl = gf.create_soil()
    sl.g_mod = 30e6
    sl.poissons_ratio = 0.3
    fd = gf.create_foundation(length, width, depth)
    i_ll = width ** 3 * length / 12
    i_ww = length ** 3 * width / 12
    k_rot = gf.stiffness.calc_rot_via_gazetas_1991(sl, fd, ip_axis='width')
    k_strip = gf.stiffness.calc_rot_strip_via_gazetas_1991(sl, fd, ip_axis='width')
    print(k_strip / k_rot)
    k_vert = gf.stiffness.calc_vert_via_gazetas_1991(sl, fd)
    k_strip = gf.stiffness.calc_vert_strip_via_gazetas_1991(sl, fd, ip_axis='width')
    print(k_strip / k_vert)
    k_horz = gf.stiffness.calc_horz_via_gazetas_1991(sl, fd, ip_axis='width')
    k_strip = gf.stiffness.calc_horz_strip_via_gazetas_1991(sl, fd, ip_axis='width')
    print(k_strip / k_horz)
