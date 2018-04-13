
import geofound


def test_schmertmann():
    """
    Checks settlement value against that from HW#7 crib ENCN452 course 2013
    """
    dry_unit_dry_weight = 18.0
    length = 2.5
    width = 2.5
    depth = 0.8
    phi = 0.0
    cohesion = 20
    unit_dry_weight = 17
    sl = geofound.create_soil(phi, cohesion, unit_dry_weight)
    fd = geofound.create_foundation(length, width, depth)
    load = 1000 + (23.5 * 0.5 + dry_unit_dry_weight * 0.3) * length * width  # kN
    youngs_modulus_soil = 3 * 6e3  # kPa
    sett = geofound.settlement_schmertmann(sl, fd, load, youngs_modulus_soil, unit_sat_weight=19.8, gwl=0.8)

    assert geofound.isclose(sett, 0.017, rel_tol=0.01)


def test_schmertmann_long_foundation():
    """
    Checks settlement value against that from xxx
    """
    dry_unit_dry_weight = 18.0
    length = 26
    width = 2.5
    depth = 0.8
    phi = 0.0
    cohesion = 20
    unit_dry_weight = 17
    sl = geofound.create_soil(phi, cohesion, unit_dry_weight)
    fd = geofound.create_foundation(length, width, depth)
    load = 1000 + (23.5 * 0.5 + dry_unit_dry_weight * 0.3) * length * width  # kN
    youngs_modulus_soil = 3 * 6e3  # kPa
    sett = geofound.settlement_schmertmann(sl, fd, load, youngs_modulus_soil, unit_sat_weight=19.8, gwl=0.8)
    print(sett)
    assert geofound.isclose(sett, 0.00153, rel_tol=0.01)


if __name__ == '__main__':
    test_schmertmann_long_foundation()