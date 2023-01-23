# -*- coding: utf-8 -*-
#  Copyright (c) 2022 Kumagai group.
from vise.tests.helpers.assertion import assert_json_roundtrip



def test_e_p_coupling_to_json_file(e_p_coupling, tmpdir):
    assert_json_roundtrip(e_p_coupling, tmpdir)


def test_reset_inner_products(e_p_coupling):
    e_p_coupling.reset_inner_prod()
    actual = e_p_coupling.e_p_matrix_elements[1][0].inner_products
    expected = []
    assert actual == expected


# def test_make_initial_e_p_coupling(dephon_init, ccd):
#     unitcell = Unitcell("", 1.0, 2.0,
#                         ele_dielectric_const=[[1.0, 0.0, 0.0],
#                                               [0.0, 2.0, 0.0],
#                                               [0.0, 0.0, 3.0]],
#                         ion_dielectric_const=[[0.0]*3]*3)
#
#     effective_mass = EffectiveMass(p=[[[10.0, 0.0, 0.0],
#                                        [0.0, 20.0, 0.0],
#                                        [0.0, 0.0, 30.0]]],
#                                    n=[[[0.0]*3]*3],
#                                    temperature=300,
#                                    concentrations=[10**18])
#
#     actual = make_initial_e_p_coupling(dephon_init=dephon_init,
#                                        ccd=ccd,
#                                        ccd_initial_name="ground + h + e",
#                                        ccd_final_name="excited + e",
#                                        unitcell=unitcell,
#                                        effective_mass=effective_mass,
#                                        base_state=BaseState.initial)
#
#     expected = EPCoupling(initial_charge=-1,
#                           initial_name="ground + h + e",
#                           final_name="excited + e",
#                           base_state=BaseState.initial,
#                           captured_carrier=Carrier.hole,
#                           volume=8.0,
#                           carrier_effective_mass=20.0,
#                           static_dielectric_const=2.0,
#                           e_p_matrix_elements=
#                           [EPMatrixElements(100,
#                                             spin_channel=Spin.up,
#                                             matrices=[EPMatrixElement(101, 1000.0, kpt_coord=[0.0, 0.0, 0.0])])])
#     assert actual == expected


"""
TODO:

"""


# def test_add_inner_prod(e_p_matrix_elements):
#     e_p_matrix_elements.add_inner_prod(band_edge_index=2,
#                                        inner_prod=2.0,
#                                        Q=20.0)
#     actual = e_p_matrix_elements.matrix_elements[0]
#     expected = EPMatrixElement(band_edge_index=2,
#                                eigenvalue_diff=0.1,
#                                inner_products=[1.0, 2.0],
#                                Qs=[10.0, 20.0],
#                                fitting_range=[-1.0, 10.1])
#     assert actual == expected


# def test_make_e_p_matrix_element(test_files):
#     dir_ = test_files / "NaP" / "Va_P1_-1_0" / "from_0_to_-1_after_make_single_point_infos"
#     actual = make_e_p_matrix_elements(
#         no_disp_dir=dir_ / "disp_0.0",
#         calc_dirs=[dir_ / "disp_0.0", dir_ / "disp_0.1"],
#         defect_index=766,
#         band_edge_indices=[767],
#         spin=Spin.down,
#         kpoint_idx=1)
#     print(actual)
# #    expected = EPMatrixElements(
#    )

