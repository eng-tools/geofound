from geofound import models as gm
from geofound import stiffness

from tests import checking_tools as ct


def test_rotation():
    length = 3
    width = 2
    depth = 1
    sl = gm.create_soil()
    sl.g_mod = 30e6
    sl.poissons_ratio = 0.3
    fd = gm.create_foundation(length, width, depth)
    i_ll = width ** 3 * length / 12
    i_ww = length ** 3 * width / 12
    assert ct.isclose(fd.i_ll, i_ll, rel_tol=0.001)
    assert ct.isclose(fd.i_ww, i_ww, rel_tol=0.001)
    assert ct.isclose(stiffness.rotational_stiffness(sl, fd, axis="length"), 229789751.26, rel_tol=0.01)
    assert ct.isclose(stiffness.rotational_stiffness(sl, fd, axis="width"), 373800747.83, rel_tol=0.01)


if __name__ == '__main__':
    test_rotation()
