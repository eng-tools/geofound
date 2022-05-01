import geofound as gf


def show_example():

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
    c_rot = gf.damping.calc_rot_via_gazetas_1991(sl, fd, ip_axis='width', a0=a0)
    c_strip = gf.damping.calc_rot_strip_via_gazetas_1991(sl, fd, ip_axis='width', a0=a0)
    print(c_strip / c_rot * fd.length)
    c_vert = gf.damping.calc_vert_via_gazetas_1991(sl, fd, a0=a0)
    c_strip = gf.damping.calc_vert_strip_via_gazetas_1991(sl, fd, ip_axis='width', a0=a0)
    print(c_strip / c_vert * fd.length)
    c_horz = gf.damping.calc_horz_via_gazetas_1991(sl, fd, ip_axis='width', a0=a0)
    c_strip = gf.damping.calc_horz_strip_via_gazetas_1991(sl, fd, ip_axis='width', a0=a0)
    print(c_strip / c_horz * fd.length)
    gf.damping.calc_rot_via_gazetas_1991(sl, fd, ip_axis='width', a0=0.5)


if __name__ == '__main__':
    show_example()
