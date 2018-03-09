from geofound import models as gm
from geofound import capacity

from geofound import checking_tools as ct


def test_vesics():
    """
    values from: Gunaratne, Manjriker. 2006. "Spread Footings: Analysis and Design."
    """
    length = 2
    width = 2
    depth = 1
    phi = 35
    cohesion = 0
    unit_dry_weight = 17
    sl = gm.create_soil(phi, cohesion, unit_dry_weight)
    fd = gm.create_foundation(length, width, depth)
    capacity.vesics_capacity(sl, fd, verbose=0)
    assert ct.isclose(fd.nc_factor, 46.1, rel_tol=0.001)
    assert ct.isclose(fd.nq_factor, 33.3, rel_tol=0.001)
    assert ct.isclose(fd.ng_factor, 48.0, rel_tol=0.001)
    assert ct.isclose(fd.q_ult, 1574.8, rel_tol=0.001)


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
    sl = gm.create_soil(phi, cohesion, unit_dry_weight)
    fd = gm.create_foundation(length, width, depth)
    capacity.terzaghi_capacity(sl, fd, verbose=0)
    assert ct.isclose(fd.nc_factor, 57.8, rel_tol=0.001)
    assert ct.isclose(fd.nq_factor, 41.4, rel_tol=0.001)
    assert ct.isclose(fd.ng_factor, 47.3, rel_tol=0.001)
    assert ct.isclose(fd.q_ult, 1347.0, rel_tol=0.001)
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
    sl = gm.create_soil(phi, cohesion, unit_dry_weight)
    fd = gm.create_foundation(length, width, depth)
    capacity.terzaghi_capacity(sl, fd, verbose=0)
    print(fd.q_ult)
    assert ct.isclose(fd.nc_factor, 12.86, rel_tol=0.001)
    assert ct.isclose(fd.nq_factor, 4.45, rel_tol=0.001)
    assert ct.isclose(fd.ng_factor, 2.168, rel_tol=0.001)
    assert ct.isclose(fd.q_ult, 515.0, rel_tol=0.001)


def test_meyerhoff():

    length = 2
    width = 2
    depth = 2
    phi = 15
    cohesion = 20.0  # kPa
    unit_dry_weight = 17  # kN/m3
    sl = gm.create_soil(phi, cohesion, unit_dry_weight)
    fd = gm.create_foundation(length, width, depth)
    capacity.meyerhoff_capacity(sl, fd, verbose=0)
    print(fd.ng_factor)
    assert ct.isclose(fd.nc_factor, 10.97, rel_tol=0.001)
    assert ct.isclose(fd.nq_factor, 3.94, rel_tol=0.01)
    assert ct.isclose(fd.ng_factor, 1.13, rel_tol=0.01)
    assert ct.isclose(fd.q_ult, 573.3, rel_tol=0.001)


def test_hansen():
    length = 2
    width = 2
    depth = 2
    phi = 15
    cohesion = 20.0  # kPa
    unit_dry_weight = 17  # kN/m3
    sl = gm.create_soil(phi, cohesion, unit_dry_weight)
    fd = gm.create_foundation(length, width, depth)
    capacity.hansen_capacity(sl, fd, verbose=0)
    print(fd.q_ult)
    assert ct.isclose(fd.nc_factor, 10.97, rel_tol=0.001)
    assert ct.isclose(fd.nq_factor, 3.94, rel_tol=0.01)
    assert ct.isclose(fd.ng_factor, 1.18, rel_tol=0.01)
    assert ct.isclose(fd.q_ult, 648.0, rel_tol=0.001)


def test_nzs_vm4():
    """
    values from: NZ Building code Clause B1 VM4 example in Appendix C
    -Retain wall example load case 1 (page 71)
    """
    length = 10000  # Actually should be a strip
    width = 2.65
    depth = 0.4
    phi = 0
    cohesion = 75.0  # kPa
    unit_dry_weight = 18  # kN/m3
    sl = gm.create_soil(phi, cohesion, unit_dry_weight)
    fd = gm.create_foundation(length, width, depth)
    h_b = 70.45 * length
    vertical_load = 131.29 * length
    h_eff_b = 1.44
    loc_v_b = 0.848
    capacity.nzs_vm4_capacity(sl, fd, h_b=h_b, vertical_load=vertical_load, h_eff_b=h_eff_b, loc_v_b=loc_v_b, verbose=0)
    print(fd.q_ult)
    assert ct.isclose(fd.nc_factor, 5.14, rel_tol=0.001)
    assert ct.isclose(fd.nq_factor, 1.0, rel_tol=0.01)
    assert ct.isclose(fd.ng_factor, 0.0, rel_tol=0.01)
    assert ct.isclose(fd.q_ult, 368.12, rel_tol=0.001)


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
    sl = gm.create_soil(phi, cohesion, unit_dry_weight)
    fd = gm.create_foundation(length, width, depth)
    h_b = 93.27 * length
    vertical_load = 154.87 * length
    h_eff_b = 1.78
    loc_v_b = 0.854
    capacity.nzs_vm4_capacity(sl, fd, h_b=h_b, vertical_load=vertical_load, h_eff_b=h_eff_b, loc_v_b=loc_v_b, verbose=0)
    print(fd.q_ult)
    assert ct.isclose(fd.nc_factor, 5.14, rel_tol=0.001)
    assert ct.isclose(fd.nq_factor, 1.0, rel_tol=0.01)
    assert ct.isclose(fd.ng_factor, 0.0, rel_tol=0.01)
    assert ct.isclose(fd.q_ult, 301.68, rel_tol=0.001)


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
    sl = gm.create_soil(phi, cohesion, unit_dry_weight)
    fd = gm.create_foundation(length, width, depth)
    h_b = 70.45 * length
    vertical_load = 144.48 * length
    h_eff_b = 1.44
    loc_v_b = 0.813
    capacity.nzs_vm4_capacity(sl, fd, h_b=h_b, vertical_load=vertical_load, h_eff_b=h_eff_b, loc_v_b=loc_v_b, verbose=0)

    assert ct.isclose(fd.nc_factor, 20.72, rel_tol=0.001)
    assert ct.isclose(fd.nq_factor, 10.66, rel_tol=0.01)
    assert ct.isclose(fd.ng_factor, 9.01, rel_tol=0.01)
    assert ct.isclose(fd.q_ult, 145.02, rel_tol=0.001)  # 152.70?


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
    sl = gm.create_soil(phi, cohesion, unit_dry_weight)
    fd = gm.create_foundation(length, width, depth)
    capacity.terzaghi_capacity(sl, fd, verbose=0)
    print(fd.nc_factor)
    print(fd.q_ult)
    # assert ct.isclose(fd.nc_factor, 6.28, rel_tol=0.001)
    # assert ct.isclose(fd.q_ult, 255.0, rel_tol=0.001)
    capacity.vesics_capacity(sl, fd, verbose=0)
    assert ct.isclose(fd.q_ult, 298.0, rel_tol=0.001)


def test_size_foundations():
    phi = 32
    cohesion = 0
    unit_dry_weight = 20
    fos_values = [3, 5, 10, 20]
    vertical_loads = [500., 800., 1000.]
    methods = capacity.available_methods
    for fos in fos_values:
        for vertical_load in vertical_loads:
            for method in methods:
                sl = gm.create_soil(phi, cohesion, unit_dry_weight)
                fd = capacity.size_footing(sl, vertical_load=vertical_load, fos=fos, method=method, length_to_width=2)
                f_capacity = fd.length * fd.width * fd.q_ult
                actual_fos = f_capacity / vertical_load
                assert ct.isclose(actual_fos, fos, rel_tol=0.1), (fos, vertical_load, method)


def test_meyerhoff_and_hanna_capacity_strong_sand_over_weak_clay():
    # STRONG SAND OVER WEAK CLAY
    length = 1000000.0  # actually a strip in
    width = 1.0
    depth = 0.0
    fd = gm.create_foundation(length=length, width=width, depth=depth)

    phi_0 = 34.0
    cohesion_0 = 0.0
    unit_dry_weight_0 = 17.0
    sl_0 = gm.create_soil(phi=phi_0, cohesion=cohesion_0, unit_dry_weight=unit_dry_weight_0)

    phi_1 = 0.0
    cohesion_1 = 30.0
    unit_dry_weight_1 = 17.0
    sl_1 = gm.create_soil(phi=phi_1, cohesion=cohesion_1, unit_dry_weight=unit_dry_weight_1)
    h0 = 1.5  # m, height of the crust layer
    #c_a = 0.0
    #k =4.8

    capacity.meyerhoff_and_hanna_capacity(sl_0, sl_1, h0, fd, verbose=0)
    # print("q_b= " + str(fd.q_b))
    print("q_ult_1= " + str(fd.q_ult))
    # print("q_t= " + str(fd.q_t))

    assert ct.isclose(fd.q_ult, 264.74, rel_tol=0.001), fd.q_ult


def test_meyerhoff_and_hanna_capacity_strong_sand_over_weak_sand():
    # STRONG SAND OVER WEAK SAND
    length = 1000000.0  # actually a strip in
    width = 1.0
    depth = 0.0
    fd = gm.create_foundation(length=length, width=width, depth=depth)

    phi_0 = 34.0
    cohesion_0 = 0.0
    unit_dry_weight_0 = 17.0
    sl_0 = gm.create_soil(phi=phi_0, cohesion=cohesion_0, unit_dry_weight=unit_dry_weight_0)

    phi_1 = 17.0
    cohesion_1 = 0.0
    unit_dry_weight_1 = 17.0
    sl_1 = gm.create_soil(phi=phi_1, cohesion=cohesion_1, unit_dry_weight=unit_dry_weight_1)
    h0 = 1.5  # m, height of the crust layer
    #c_a = 0.0
    #k=1.7

    capacity.meyerhoff_and_hanna_capacity(sl_0, sl_1, h0, fd, verbose=0)
    #print("q_b= " +str(fd.q_b))
    #print("q_ult5= " + str(fd.q_ult5))
    print("q_ult_2= " +str(fd.q_ult))
    assert ct.isclose(fd.q_ult, 158.32, rel_tol=0.001), fd.q_ult


def test_meyerhoff_and_hanna_capacity_strong_clay_over_weak_sand():
    # STRONG CLAY OVER WEAK SAND
    length = 1000000000.0  # actually a strip in
    width = 1.0
    depth = 0.0
    fd = gm.create_foundation(length=length, width=width, depth=depth)

    phi_0 = 0.0
    cohesion_0 = 85.0
    unit_dry_weight_0 = 16.5
    sl_0 = gm.create_soil(phi=phi_0, cohesion=cohesion_0, unit_dry_weight=unit_dry_weight_0)

    phi_1 = 17.0
    cohesion_1 = 0.0
    unit_dry_weight_1 = 17.0
    sl_1 = gm.create_soil(phi=phi_1, cohesion=cohesion_1, unit_dry_weight=unit_dry_weight_1)
    h0 = 1  # m, height of the crust layer
    #c_a = 54.4
    #k = 1.7

    capacity.meyerhoff_and_hanna_capacity(sl_0, sl_1, h0, fd, verbose=0)
    #print(sl_0.nc_factor_0)
    #print(fd.q_ult)

    #print("q_b= " + str(fd.q_b))
    #print("q1_q0= " + str(fd.q1_q0))
    print("q_ult_3= " + str(fd.q_ult))
    assert ct.isclose(fd.q_ult, 182.97, rel_tol=0.001), fd.q_ult

# def test_meyerhoff_and_hanna_capacity_strong_clay_over_weak_sand_vs_limitstategeo():
#     # STRONG CLAY OVER WEAK SAND
#     length = 1000000.0  # actually a strip in
#     width = 10.0
#     depth = 0.0
#     fd = gm.create_foundation(length=length, width=width, depth=depth)
#
#     phi_0 = 0.0
#     cohesion_0 = 50.0
#     unit_dry_weight_0 = 19.
#     sl_0 = gm.create_soil(phi=phi_0, cohesion=cohesion_0, unit_dry_weight=unit_dry_weight_0)
#
#     phi_1 = 3.
#     cohesion_1 = 0.0
#     unit_dry_weight_1 = 19.
#     sl_1 = gm.create_soil(phi=phi_1, cohesion=cohesion_1, unit_dry_weight=unit_dry_weight_1)
#     h0 = 2.  # m, height of the crust layer
#     #c_a = 54.4
#     #k = 1.7
#
#     capacity.meyerhoff_and_hanna_capacity(sl_0, sl_1, h0, fd, verbose=0)
#     #print(sl_0.nc_factor_0)
#     #print(fd.q_ult)
#
#     #print("q_b= " + str(fd.q_b))
#     print("q_b_limitstategeo= " + str(fd.q_b))
#     print("q_t_limitstategeo= " + str(fd.q_t))
#     print("q_ult_limitstategeo= " + str(fd.q_ult))
#     assert ct.isclose(fd.q_ult, 26.45, rel_tol=0.001), fd.q_ult

def load_soil_sample_data(sl0):
    """
    Sample data for the Soil object
    :param sl0: Soil Object
    :return:
    """
    # soil
    sl0.height = 1.5 #[m]
    sl0.phi = 34  # [degrees]
    sl0.unit_dry_weight = 17000  # [N/m3]
    sl0.c_a = 0
    sl0.cohesion = 0 # [Pa]

def load_soil_sample_data2(sl1):
    """
    Sample data for the Soil object
    :param sl1: Soil Object
    :return:
    """
    # soil
    sl1.cohesion = 30 # [Pa]
    sl1.phi = 0  # [degrees]
    sl1.unit_dry_weight = 17000  # [N/m3]



if __name__ == '__main__':
    test_meyerhoff_and_hanna_capacity_strong_clay_over_weak_sand()