import geofound
import numpy as np
from geofound import models


def test_vesic():
    """
    values from: Gunaratne, Manjriker. 2006. "Spread Footings: Analysis and Design."
    """
    length = 2
    width = 2
    depth = 1
    phi = 35
    cohesion = 0
    unit_dry_weight = 17
    sl = geofound.create_soil(phi, cohesion, unit_dry_weight)
    fd = geofound.create_foundation(length, width, depth)
    geofound.capacity_vesic_1975(sl, fd, verbose=0)
    assert np.isclose(fd.nc_factor, 46.1, rtol=0.001)
    assert np.isclose(fd.nq_factor, 33.3, rtol=0.001)
    assert np.isclose(fd.ng_factor, 48.0, rtol=0.001)
    assert np.isclose(fd.q_ult, 1574.8, rtol=0.001)


def test_terzaghi():
    """
    values from: Gunaratne, Manjriker. 2006. "Spread Footings: Analysis and Design."
    - except qult was not validated
    """
    length = 2
    width = 2
    depth = 1
    phi = 35
    cohesion = 0
    unit_dry_weight = 17
    sl = geofound.create_soil(phi, cohesion, unit_dry_weight)
    fd = geofound.create_foundation(length, width, depth)
    geofound.capacity_terzaghi_1943(sl, fd, verbose=0)
    assert np.isclose(fd.nc_factor, 57.8, rtol=0.001)
    assert np.isclose(fd.nq_factor, 41.4, rtol=0.001)
    assert np.isclose(fd.ng_factor, 47.3, rtol=0.001)
    assert np.isclose(fd.q_ult, 1347.0, rtol=0.001)
    print("DONE")


def test_terzaghi_again():
    """
    values from: Gunaratne, Manjriker. 2006. "Spread Footings: Analysis and Design."
    """
    length = 2
    width = 2
    depth = 2
    phi = 15
    cohesion = 20.0  # kPa
    unit_dry_weight = 17  # kN/m3
    sl = geofound.create_soil(phi, cohesion, unit_dry_weight)
    fd = geofound.create_foundation(length, width, depth)
    geofound.capacity_terzaghi_1943(sl, fd, verbose=0)
    print(fd.q_ult)
    assert np.isclose(fd.nc_factor, 12.86, rtol=0.001)
    assert np.isclose(fd.nq_factor, 4.45, rtol=0.001)
    assert np.isclose(fd.ng_factor, 2.168, rtol=0.001)
    assert np.isclose(fd.q_ult, 515.0, rtol=0.001)


def test_meyerhof():
    length = 2
    width = 2
    depth = 2
    phi = 15
    cohesion = 20.0  # kPa
    unit_dry_weight = 17  # kN/m3
    sl = geofound.create_soil(phi, cohesion, unit_dry_weight)
    fd = geofound.create_foundation(length, width, depth)
    geofound.capacity_meyerhof_1963(sl, fd, gwl=1000, verbose=0)
    print(fd.ng_factor)
    assert np.isclose(fd.nc_factor, 10.97, rtol=0.001)
    assert np.isclose(fd.nq_factor, 3.94, rtol=0.01)
    assert np.isclose(fd.ng_factor, 1.13, rtol=0.01)
    assert np.isclose(fd.q_ult, 573.3, rtol=0.001)


def test_salgado_2008():
    length = 3
    width = 2
    depth = 0.75
    phi = 36.1
    cohesion = 0.0  # kPa
    unit_dry_weight = 20.0  # kN/m3
    sl = geofound.create_soil(phi, cohesion, unit_dry_weight)
    fd = geofound.create_foundation(length, width, depth)
    geofound.capacity_salgado_2008(sl, fd, gwl=1000, verbose=0, save_factors=1)
    print(fd.ng_factor)
    assert np.isclose(fd.nq_factor, 38.3, rtol=0.01)  # page 446
    assert np.isclose(fd.ng_factor, 40.9, rtol=0.01)
    assert np.isclose(fd.q_d, 15.0, rtol=0.001), fd.qd
    assert np.isclose(fd.s_q, 1.935, rtol=0.01), getattr(fd, 's_q')
    assert np.isclose(fd.s_g, 1.12, rtol=0.01), fd.s_g
    assert np.isclose(fd.d_q, 1.68, rtol=0.01), fd.d_q
    assert np.isclose(fd.d_g, 1.0, rtol=0.01), fd.d_g

    assert np.isclose(fd.q_ult, 2780.0, rtol=0.001), fd.q_ult  # Note: differs from text due to rounding in the text


def test_meyerhof_using_fabrizio_problem1():
    """
    http:

    :return:
    """

    length = 100000
    width = 4.
    depth = 1.5
    phi = 0.0
    cohesion = 90000  # Pa
    unit_dry_weight = 19000.  # N/m3
    unit_sat_weight = 19000.
    sl = geofound.create_soil(phi, cohesion, unit_dry_weight)
    sl.unit_sat_weight = unit_sat_weight
    fd = geofound.create_foundation(length, width, depth)
    geofound.capacity_meyerhof_1963(sl, fd, verbose=0)
    assert np.isclose(fd.q_ult, 526000., rtol=0.001), fd.q_ult


def test_meyerhof_using_fabrizio_problem2():
    """
    http:

    :return:
    """

    length = 100000
    width = 4.
    depth = 1.5
    phi = 34.0
    cohesion = 0.0  # kPa
    unit_dry_weight = 18.  # kN/m3
    unit_sat_weight = 20.
    sl = geofound.create_soil(phi, cohesion, unit_dry_weight, pw=9.8)

    sl.unit_sat_weight = unit_sat_weight
    fd = geofound.create_foundation(length, width, depth)

    # problem 2) a)

    geofound.capacity_meyerhof_1963(sl, fd, gwl=20.0, verbose=0)
    # # assert np.isclose(fd.nc_factor, 10.97, rtol=0.001)
    assert np.isclose(fd.nq_factor, 29.4, rtol=0.01)
    assert np.isclose(fd.ng_factor, 31.1, rtol=0.01)
    assert np.isclose(fd.q_ult, 2056, rtol=0.01)

    # problem 2) b)
    geofound.capacity_meyerhof_1963(sl, fd, gwl=1.5, verbose=0)
    # assert np.isclose(fd.nc_factor, 10.97, rtol=0.001)
    assert np.isclose(fd.nq_factor, 29.4, rtol=0.01)
    assert np.isclose(fd.ng_factor, 31.1, rtol=0.01)
    assert np.isclose(fd.q_ult, 1521, rtol=0.01), fd.q_ult

    # problem 2) c)
    geofound.capacity_meyerhof_1963(sl, fd, gwl=0.5, verbose=0)
    # assert np.isclose(fd.nc_factor, 10.97, rtol=0.001)
    assert np.isclose(fd.nq_factor, 29.4, rtol=0.01)
    assert np.isclose(fd.ng_factor, 31.1, rtol=0.01)
    assert np.isclose(fd.q_ult, 1252, rtol=0.03), fd.q_ult


def test_hansen():
    length = 2
    width = 2
    depth = 2
    phi = 15
    cohesion = 20.0  # kPa
    unit_dry_weight = 17  # kN/m3
    sl = geofound.create_soil(phi, cohesion, unit_dry_weight)
    fd = geofound.create_foundation(length, width, depth)
    geofound.capacity_brinch_hansen_1970(sl, fd, verbose=0)
    assert np.isclose(fd.nc_factor, 10.97, rtol=0.001)
    assert np.isclose(fd.nq_factor, 3.94, rtol=0.01)
    assert np.isclose(fd.ng_factor, 1.18, rtol=0.01)
    assert np.isclose(fd.q_ult, 648.0, rtol=0.001)


def test_nzs_vm4():
    """
    values from: NZ Building code Clause B1 VM4 example in Appendix C
    - Retain wall example load case 1 (page 71)
    """
    length = 1  # Actually should be a strip
    width = 2.65
    depth = 0.4
    phi = 0
    cohesion = 75.0  # kPa
    unit_dry_weight = 18  # kN/m3
    sl = geofound.create_soil(phi, cohesion, unit_dry_weight)
    fd = geofound.create_foundation(length, width, depth)
    h_b = 70.47 * length
    vertical_load = 154.87 * length
    h_eff_b = 1.44
    loc_v_b = 0.848  # X
    e_applied = width / 2 - loc_v_b
    mom_width = h_b * h_eff_b
    e_mom = mom_width / vertical_load
    e_width = e_mom - e_applied
    geofound.capacity_nzs_vm4_2011(sl, fd, hload_width=h_b, nload=vertical_load, e_width=e_width,
                                   ip_axis_2d='width', save_factors=1, verbose=0)
    print('width_eff: ', fd.width_eff)
    assert np.isclose(fd.width_eff, 2.29, rtol=0.01), fd.width_eff
    assert np.isclose(fd.nc_factor, 5.14, rtol=0.001)
    assert np.isclose(fd.nq_factor, 1.0, rtol=0.01)
    assert np.isclose(fd.ng_factor, 0.0, rtol=0.01)
    assert np.isclose(fd.d_c, 1.07, rtol=0.01)
    assert np.isclose(fd.s_c, 1.0, rtol=0.0001)
    assert np.isclose(fd.i_c, 0.88, rtol=0.01)
    assert np.isclose(fd.q_ult, 371.8227, rtol=0.001), fd.q_ult  # Note reported answer is 370.19 due to rounding



def test_nzs_vm4_load_case_3():
    """
    values from: NZ Building code Clause B1 VM4 example in Appendix C
    -Retain wall example load case 3
    """
    length = 10000  # Actually should be a strip
    width = 2.65
    depth = 0.4
    phi = 0
    cohesion = 75.0  # kPa
    unit_dry_weight = 18  # kN/m3
    sl = geofound.create_soil(phi, cohesion, unit_dry_weight)
    fd = geofound.create_foundation(length, width, depth)
    h_b = 93.27 * length
    vertical_load = 154.87 * length
    h_eff_b = 1.78
    loc_v_b = 0.854
    e_applied = width / 2 - loc_v_b
    mom_width = h_b * h_eff_b
    e_mom = mom_width / vertical_load
    e_width = e_mom - e_applied
    geofound.capacity_nzs_vm4_2011(sl, fd, hload_width=h_b, nload=vertical_load, e_width=e_width,
                                   verbose=0)
    assert np.isclose(fd.nc_factor, 5.14, rtol=0.001)
    assert np.isclose(fd.nq_factor, 1.0, rtol=0.01)
    assert np.isclose(fd.ng_factor, 0.0, rtol=0.01)
    assert np.isclose(fd.q_ult, 301.68, rtol=0.001)


def test_nzs_vm4_load_case_5():
    """
    values from: NZ Building code Clause B1 VM4 example in Appendix C
    -Retain wall example load case 5
    ***looks like there is an error in the calculation of d_c
    ***small discrepancy between Nc, i_c and i_q
    """
    length = 10000  # Actually should be a strip
    width = 2.65
    depth = 0.4
    phi = 25.0
    cohesion = 12.5  # kPa
    unit_dry_weight = 8.2  # kN/m3
    sl = geofound.create_soil(phi, cohesion, unit_dry_weight)
    fd = geofound.create_foundation(length, width, depth)
    h_b = 70.45 * length
    vertical_load = 144.48 * length
    h_eff_b = 1.44
    loc_v_b = 0.813
    e_applied = width / 2 - loc_v_b
    mom_width = h_b * h_eff_b
    e_mom = mom_width / vertical_load
    e_width = e_mom - e_applied
    geofound.capacity_nzs_vm4_2011(sl, fd, hload_width=h_b, nload=vertical_load, e_width=e_width,
                                   verbose=0)
    assert np.isclose(fd.nc_factor, 20.72, rtol=0.001)
    assert np.isclose(fd.nq_factor, 10.66, rtol=0.01)
    assert np.isclose(fd.ng_factor, 9.01, rtol=0.01)
    assert np.isclose(fd.q_ult, 145.02, rtol=0.001)  # 152.70?


def test_from_encn452_2013():
    """
    Values from HW#7 crib ENCN452 course 2013
    """
    length = 6.0  # actually a strip in
    width = 3.0
    depth = 1.5
    phi = 0.0
    cohesion = 40.0
    unit_dry_weight = 18.0
    sl = geofound.create_soil(phi, cohesion, unit_dry_weight)
    fd = geofound.create_foundation(length, width, depth)
    geofound.capacity_terzaghi_1943(sl, fd, verbose=0)
    geofound.capacity_vesic_1975(sl, fd, verbose=0)
    assert np.isclose(fd.q_ult, 298.0, rtol=0.001)


def test_size_foundations():
    phi = 32
    cohesion = 0
    unit_dry_weight = 20
    fos_values = [3.0, 5.0, 10.0, 20.0]
    vertical_loads = [500., 800., 1000.]
    methods = geofound.available_methods
    for fos in fos_values:
        for vertical_load in vertical_loads:
            for method in methods:
                sl = geofound.create_soil(phi, cohesion, unit_dry_weight)
                fd = geofound.size_footing_for_capacity(sl, vertical_load=vertical_load, fos=fos, method=method,
                                                        length_to_width=2)
                f_capacity = fd.length * fd.width * fd.q_ult
                actual_fos = f_capacity / vertical_load
                assert np.isclose(actual_fos, fos, rtol=0.11), (fos, vertical_load, method)


def test_calc_crit_length():
    length = 6.0  # actually a strip in
    width = 3.0
    depth = 1.5
    phi = 0.0
    cohesion = 40.0
    unit_dry_weight = 18.0
    sl = geofound.create_soil(phi, cohesion, unit_dry_weight)
    fd = geofound.create_foundation(length, width, depth)
    # vload = 3000.0
    q_ult = geofound.capacity_vesic_1975(sl, fd)
    vload = 0.5 * q_ult * fd.area
    crit_len = geofound.capacity.calc_crit_span(sl, fd, vload, ip_axis='length')

    assert np.isclose(crit_len, 2.765625, rtol=0.001), crit_len


def test_very_short_foundation_w_vesic():
    length = 0.265  # actually a strip in
    width = 4.7
    depth = 2.24
    fd = geofound.create_foundation(length, width, depth)
    width = 0.265
    length = 4.7
    fdo = geofound.create_foundation(length, width, depth)
    phi = 32.0
    cohesion = 0.0
    unit_dry_weight = 14943.15
    sl = geofound.create_soil(phi, cohesion, unit_dry_weight)

    q_ult = geofound.capacity_vesic_1975(sl, fd)
    q_ulto = geofound.capacity_vesic_1975(sl, fdo)
    assert np.isclose(q_ult, q_ulto, rtol=0.001)


def test_meyerhof_and_hanna_capacity_strong_sand_over_weak_clay():
    # STRONG SAND OVER WEAK CLAY
    length = 1000000.0  # actually a strip in
    width = 1.0
    depth = 0.0
    fd = geofound.create_foundation(length=length, width=width, depth=depth)

    phi_0 = 34.0
    cohesion_0 = 0.0
    unit_dry_weight_0 = 17.0
    sl_0 = geofound.create_soil(phi=phi_0, cohesion=cohesion_0, unit_dry_weight=unit_dry_weight_0)

    phi_1 = 0.0
    cohesion_1 = 30.0
    unit_dry_weight_1 = 17.0
    sl_1 = geofound.create_soil(phi=phi_1, cohesion=cohesion_1, unit_dry_weight=unit_dry_weight_1)
    h0 = 1.5  # m, height of the crust layer

    geofound.capacity_meyerhof_and_hanna_1978(sl_0, sl_1, h0, fd, verbose=0)

    assert np.isclose(fd.q_ult, 264.74, rtol=0.001), fd.q_ult


def test_meyerhof_and_hanna_capacity_strong_sand_over_weak_sand():
    # STRONG SAND OVER WEAK SAND
    length = 1000000.0  # actually a strip in
    width = 1.0
    depth = 0.0
    fd = geofound.create_foundation(length=length, width=width, depth=depth)

    phi_0 = 34.0
    cohesion_0 = 0.0
    unit_dry_weight_0 = 17.0
    sl_0 = geofound.create_soil(phi=phi_0, cohesion=cohesion_0, unit_dry_weight=unit_dry_weight_0)

    phi_1 = 17.0
    cohesion_1 = 0.0
    unit_dry_weight_1 = 17.0
    sl_1 = geofound.create_soil(phi=phi_1, cohesion=cohesion_1, unit_dry_weight=unit_dry_weight_1)
    h0 = 1.5  # m, height of the crust layer

    geofound.capacity_meyerhof_and_hanna_1978(sl_0, sl_1, h0, fd, verbose=0)
    assert np.isclose(fd.q_ult, 158.32, rtol=0.001), fd.q_ult


def test_meyerhof_and_hanna_capacity_strong_clay_over_weak_sand():
    # STRONG CLAY OVER WEAK SAND
    length = 1000000000.0  # actually a strip in
    width = 1.0
    depth = 0.0
    fd = geofound.create_foundation(length=length, width=width, depth=depth)

    phi_0 = 0.0
    cohesion_0 = 85.0
    unit_dry_weight_0 = 16.5
    sl_0 = geofound.create_soil(phi=phi_0, cohesion=cohesion_0, unit_dry_weight=unit_dry_weight_0)

    phi_1 = 17.0
    cohesion_1 = 0.0
    unit_dry_weight_1 = 17.0
    sl_1 = geofound.create_soil(phi=phi_1, cohesion=cohesion_1, unit_dry_weight=unit_dry_weight_1)
    h0 = 1  # m, height of the crust layer

    geofound.capacity_meyerhof_and_hanna_1978(sl_0, sl_1, h0, fd, verbose=0)
    assert np.isclose(fd.q_ult, 187.87, rtol=0.001), fd.q_ult  # no independent validation


def test_meyerhof_and_hanna_capacity_sand_over_sand_gwl():
    # STRONG SAND OVER WEAK SAND
    length = 1000000.0  # actually a strip in
    width = 4.0
    depth = 1.5
    fd = geofound.create_foundation(length=length, width=width, depth=depth)

    phi_0 = 34.0
    cohesion_0 = 0.0
    unit_dry_weight_0 = 18000
    sl_0 = geofound.create_soil(phi=phi_0, cohesion=cohesion_0, unit_dry_weight=unit_dry_weight_0)
    sl_0.unit_sat_weight = 20000

    phi_1 = 34.0
    cohesion_1 = 0.0
    unit_dry_weight_1 = 18000
    sl_1 = geofound.create_soil(phi=phi_1, cohesion=cohesion_1, unit_dry_weight=unit_dry_weight_1)
    sl_1.unit_sat_weight = 20000
    h0 = 3.0  # m, height of the crust layer

    # Case 1: GWL at surface
    gwl = 0.0
    q_ult_meyerhof = geofound.capacity_meyerhof_1963(sl_0, fd, gwl=gwl)
    geofound.capacity_meyerhof_and_hanna_1978(sl_0, sl_1, h0, fd, gwl=gwl, verbose=0)
    correction = 1.07
    corrected_2layer = fd.q_ult * correction
    assert np.isclose(corrected_2layer, q_ult_meyerhof, rtol=0.01), (corrected_2layer, q_ult_meyerhof / 1000)

    # Case 2: GWL at between foundation depth and surface
    gwl = 0.5
    q_ult_meyerhof = geofound.capacity_meyerhof_1963(sl_0, fd, gwl=gwl)
    geofound.capacity_meyerhof_and_hanna_1978(sl_0, sl_1, h0, fd, gwl=gwl, verbose=0)
    correction = 1.07
    corrected_2layer = fd.q_ult * correction
    assert np.isclose(corrected_2layer, q_ult_meyerhof, rtol=0.01), (corrected_2layer, q_ult_meyerhof / 1000)

    # Case 3: GWL at between foundation depth and foundation depth plus width, and GWL < layer 1 depth
    gwl = 1.8
    q_ult_meyerhof = geofound.capacity_meyerhof_1963(sl_0, fd, gwl=gwl)
    geofound.capacity_meyerhof_and_hanna_1978(sl_0, sl_1, h0, fd, gwl=gwl, verbose=0)
    correction = 1.07
    corrected_2layer = fd.q_ult * correction
    assert np.isclose(corrected_2layer, q_ult_meyerhof, rtol=0.01), (corrected_2layer, q_ult_meyerhof / 1000)

    # Case 4: GWL at between foundation depth and foundation depth plus width, and GWL > layer 1 depth
    gwl = 4.8
    q_ult_meyerhof = geofound.capacity_meyerhof_1963(sl_0, fd, gwl=gwl)
    geofound.capacity_meyerhof_and_hanna_1978(sl_0, sl_1, h0, fd, gwl=gwl, verbose=0)
    correction = 1.03
    corrected_2layer = fd.q_ult * correction
    assert np.isclose(corrected_2layer, q_ult_meyerhof, rtol=0.01), (corrected_2layer, q_ult_meyerhof / 1000)

    # Case 5: GWL beyond foundation depth plus width
    gwl = 20.
    q_ult_meyerhof = geofound.capacity_meyerhof_1963(sl_0, fd, gwl=gwl)
    geofound.capacity_meyerhof_and_hanna_1978(sl_0, sl_1, h0, fd, gwl=gwl, verbose=0)
    correction = 1.07
    corrected_2layer = fd.q_ult * correction
    assert np.isclose(corrected_2layer, q_ult_meyerhof, rtol=0.01), (corrected_2layer, q_ult_meyerhof/1000)


def test_meyerhof_and_hanna_capacity_clay_over_clay_gwl():
    length = 1000000.0  # actually a strip in
    width = 4.0
    depth = 1.5
    fd = geofound.create_foundation(length=length, width=width, depth=depth)

    phi_0 = 0.0
    cohesion_0 = 40.0
    unit_dry_weight_0 = 18000
    sl_0 = geofound.create_soil(phi=phi_0, cohesion=cohesion_0, unit_dry_weight=unit_dry_weight_0)
    sl_0.unit_sat_weight = 20000

    phi_1 = 0.0
    cohesion_1 = 40.0
    unit_dry_weight_1 = 18000
    sl_1 = geofound.create_soil(phi=phi_1, cohesion=cohesion_1, unit_dry_weight=unit_dry_weight_1)
    sl_1.unit_sat_weight = 20000
    h0 = 3.0  # m, height of the crust layer

    # Case 1: GWL at surface
    gwl = 0.0
    q_ult_meyerhof = geofound.capacity_meyerhof_1963(sl_0, fd, gwl=gwl)
    geofound.capacity_meyerhof_and_hanna_1978(sl_0, sl_1, h0, fd, gwl=gwl, verbose=0)
    correction = 1.00
    corrected_2layer = fd.q_ult * correction
    assert np.isclose(corrected_2layer, q_ult_meyerhof, rtol=0.01), (corrected_2layer, q_ult_meyerhof / 1000)

    # Case 2: GWL at between foundation depth and surface
    gwl = 0.5
    q_ult_meyerhof = geofound.capacity_meyerhof_1963(sl_0, fd, gwl=gwl)
    geofound.capacity_meyerhof_and_hanna_1978(sl_0, sl_1, h0, fd, gwl=gwl, verbose=0)
    correction = 1.0
    corrected_2layer = fd.q_ult * correction
    assert np.isclose(corrected_2layer, q_ult_meyerhof, rtol=0.01), (corrected_2layer, q_ult_meyerhof / 1000)

    # Case 3: GWL at between foundation depth and foundation depth plus width, and GWL < layer 1 depth
    gwl = 1.8
    q_ult_meyerhof = geofound.capacity_meyerhof_1963(sl_0, fd, gwl=gwl)
    geofound.capacity_meyerhof_and_hanna_1978(sl_0, sl_1, h0, fd, gwl=gwl, verbose=0)
    correction = 1.0
    corrected_2layer = fd.q_ult * correction
    assert np.isclose(corrected_2layer, q_ult_meyerhof, rtol=0.01), (corrected_2layer, q_ult_meyerhof / 1000)

    # Case 4: GWL at between foundation depth and foundation depth plus width, and GWL > layer 1 depth
    gwl = 4.8
    q_ult_meyerhof = geofound.capacity_meyerhof_1963(sl_0, fd, gwl=gwl)
    geofound.capacity_meyerhof_and_hanna_1978(sl_0, sl_1, h0, fd, gwl=gwl, verbose=0)
    correction = 1.0
    corrected_2layer = fd.q_ult * correction
    assert np.isclose(corrected_2layer, q_ult_meyerhof, rtol=0.01), (corrected_2layer, q_ult_meyerhof / 1000)

    # Case 5: GWL beyond foundation depth plus width
    gwl = 20.
    q_ult_meyerhof = geofound.capacity_meyerhof_1963(sl_0, fd, gwl=gwl)
    geofound.capacity_meyerhof_and_hanna_1978(sl_0, sl_1, h0, fd, gwl=gwl, verbose=0)
    correction = 1.0
    corrected_2layer = fd.q_ult * correction
    assert np.isclose(corrected_2layer, q_ult_meyerhof, rtol=0.01), (corrected_2layer, q_ult_meyerhof/1000)


def test_capacity_sp_meyerhof_and_hanna_1978():
    length = 1000000.0  # actually a strip in
    width = 16.0
    depth = 0.0
    fd = geofound.create_foundation(length=length, width=width, depth=depth)

    phi_0 = 0.0
    cohesion_0 = 50000
    unit_dry_weight_0 = 18000
    sl_0 = geofound.create_soil(phi=phi_0, cohesion=cohesion_0, unit_dry_weight=unit_dry_weight_0)
    sl_0.unit_sat_weight = 20000

    phi_1 = 25
    cohesion_1 = 0.0
    unit_dry_weight_1 = 18000
    sl_1 = geofound.create_soil(phi=phi_1, cohesion=cohesion_1, unit_dry_weight=unit_dry_weight_1)
    sl_1.unit_sat_weight = 20000
    h0 = 2.0  # m, height of the crust layer
    # c_a = 0.0
    # k=1.7

    # Case 1: GWL at surface
    gwl = 100.0
    sp = models.SoilProfile()
    sp.add_layer(0, sl_0)
    sp.add_layer(h0, sl_1)
    sp.gwl = gwl
    q_ult_meyerhof = geofound.capacity_meyerhof_1963(sl_0, fd, gwl=gwl)
    # q_ult_meyerhof_and_hanna = geofound.capacity_meyerhof_and_hanna_1978(sl_0, sl_1, h0, fd, gwl=gwl, verbose=0)
    geofound.capacity_sp_meyerhof_and_hanna_1978(sp, fd, verbose=0)
    correction_lower_layer = 1.1  # unconfirmed value, test added for v0.4.6
    expected = q_ult_meyerhof * correction_lower_layer
    assert np.isclose(expected, fd.q_ult, rtol=0.01), (expected / 1000, fd.q_ult / 1000)


def test_calc_m_eff_via_loukidis_and_salgado_2006():
    # pg 451
    fd = models.RaftFoundation()
    fd.width = 2
    fd.length = 3
    fd.depth = 0.0
    sl = models.Soil(unit_dry_weight=17.1e3)
    m_eff_expected = 742.0  # kPa
    m_eff = geofound.calc_m_eff_bc_via_loukidis_and_salgado_2006(sl, fd, p_atm=100.0e3) / 1e3
    assert np.isclose(m_eff, m_eff_expected, rtol=0.01), m_eff

# def test_meyerhof_and_hanna_capacity_strong_clay_over_weak_sand_vs_limitstategeo():
# # STRONG CLAY OVER WEAK SAND
# length = 1000000.0  # actually a strip in
#     width = 10.0
#     depth = 0.0
#     fd = geofound.create_foundation(length=length, width=width, depth=depth)
#
#     phi_0 = 0.0
#     cohesion_0 = 50.0
#     unit_dry_weight_0 = 19.
#     sl_0 = geofound.create_soil(phi=phi_0, cohesion=cohesion_0, unit_dry_weight=unit_dry_weight_0)
#
#     phi_1 = 3.
#     cohesion_1 = 0.0
#     unit_dry_weight_1 = 19.
#     sl_1 = geofound.create_soil(phi=phi_1, cohesion=cohesion_1, unit_dry_weight=unit_dry_weight_1)
#     h0 = 2.  # m, height of the crust layer
#     #c_a = 54.4
#     #k = 1.7
#
#     geofound.capacity_meyerhof_and_hanna_1978(sl_0, sl_1, h0, fd, verbose=0)
#     #print(sl_0.nc_factor_0)
#     #print(fd.q_ult)
#
#     #print("q_b= " + str(fd.q_b))
#     print("q_b_limitstategeo= " + str(fd.q_b))
#     print("q_t_limitstategeo= " + str(fd.q_t))
#     print("q_ult_limitstategeo= " + str(fd.q_ult))
#     assert np.isclose(fd.q_ult, 26.45, rtol=0.001), fd.q_ult

def load_soil_sample_data(sl0):
    """
    Sample data for the Soil object
    :param sl0: Soil Object
    :return:
    """
    # soil
    sl0.height = 1.5  # [m]
    sl0.phi = 34  # [degrees]
    sl0.unit_dry_weight = 17000  # [N/m3]
    sl0.c_a = 0
    sl0.cohesion = 0  # [Pa]


def load_soil_sample_data2(sl1):
    """
    Sample data for the Soil object
    :param sl1: Soil Object
    :return:
    """
    # soil
    sl1.cohesion = 30  # [Pa]
    sl1.phi = 0  # [degrees]
    sl1.unit_dry_weight = 17000  # [N/m3]


if __name__ == '__main__':
    test_calc_crit_length()
    # test_meyerhof_and_hanna_capacity_sand_over_sand_gwl()
    # test_meyerhof_and_hanna_capacity_strong_clay_over_weak_sand()
