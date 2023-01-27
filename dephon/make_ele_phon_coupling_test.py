# -*- coding: utf-8 -*-
#  Copyright (c) 2022 Kumagai group.
from copy import deepcopy

import pytest
from pydefect.analyzer.band_edge_states import LocalizedOrbital
from pymatgen.core import Lattice, Structure

from dephon.config_coord import Ccd, SingleCcd, SinglePointInfo, SingleCcdId
from dephon.dephon_init import DephonInit, MinimumPointInfo, NearEdgeState
from dephon.ele_phon_coupling import EPCoupling, InnerProduct, EPMatrixElement
from dephon.enum import Carrier
from dephon.make_ele_phon_coupling import MakeInitialEPCoupling, \
    add_inner_products


@pytest.fixture
def dephon_init():
    lattice = Lattice.cubic(2.0)
    coords = [[0.0, 0.0, 0.0]]
    excited = MinimumPointInfo(charge=0,
                               structure=Structure(lattice, species=["H"], coords=coords),
                               energy=1.0,
                               correction_energy=1.0,
                               magnetization=-1.0,
                               localized_orbitals=[[]],
                               initial_site_symmetry="1",
                               final_site_symmetry="1",
                               parsed_dir="/",
                               valence_bands=[[NearEdgeState(band_index=101, kpt_coord=[0.0]*3, kpt_weight=1.0, kpt_index=1, eigenvalue=1000.0, occupation=0.0)], []],
                               conduction_bands=[[], []])
    ground = deepcopy(excited)
    ground.charge = -1
    return DephonInit(defect_name="test",
                      states=[ground, excited],
                      vbm=0.0, cbm=10.0,
                      supercell_vbm=1.0,
                      supercell_cbm=9.0,
                      ave_electron_mass=0.01,
                      ave_hole_mass=0.02,
                      ave_static_diele_const=0.03)


@pytest.fixture
def ccd():
    l_orb_up = LocalizedOrbital(
        band_idx=102, ave_energy=1001.0, occupation=1.0, orbitals={})
    l_orb_down = LocalizedOrbital(
        band_idx=103, ave_energy=1002.0, occupation=0.0, orbitals={})
    single_ccd_1 = SingleCcd(SingleCcdId("excited", carriers=[Carrier.e]),
                             charge=0,
                             point_infos=[SinglePointInfo(dQ=1.0, disp_ratio=0.0, magnetization=-1.0, localized_orbitals=[[l_orb_up, l_orb_down], []])])
    single_ccd_2 = SingleCcd(SingleCcdId("ground", carriers=[Carrier.e, Carrier.h]),
                             charge=-1,
                             point_infos=[SinglePointInfo(dQ=0.0, disp_ratio=0.0, magnetization=0.0, localized_orbitals=[[], []])])
    return Ccd("test", ccds=[single_ccd_1, single_ccd_2])


@pytest.fixture
def make_e_p_coupling_h(dephon_init, ccd):
    return MakeInitialEPCoupling(dephon_init=dephon_init, ccd=ccd, captured_carrier=Carrier.h)


@pytest.fixture
def make_e_p_coupling_e(dephon_init, ccd):
    return MakeInitialEPCoupling(dephon_init=dephon_init, ccd=ccd, captured_carrier=Carrier.e)


def test_make_initial_e_p_coupling_base_charge(make_e_p_coupling_h, dephon_init, ccd):
    assert make_e_p_coupling_h.charge == 0

    common = dict(dephon_init=dephon_init, ccd=ccd, captured_carrier=Carrier.h)
    make_coupling = MakeInitialEPCoupling(**common, charge_for_e_p_coupling=-1)
    assert make_coupling.charge == -1

    with pytest.raises(ValueError):
        MakeInitialEPCoupling(**common, charge_for_e_p_coupling=-2)


def test_make_initial_e_p_coupling(make_e_p_coupling_h):
    e_p_elem = EPMatrixElement(band_edge_index=101,
                               defect_band_index=102,
                               spin="up",
                               eigenvalue_diff=1.0,
                               kpt_idx=1,
                               kpt_coord=[0.0]*3)
    expected = EPCoupling(charge=0,
                          disp=0.0,
                          captured_carrier=Carrier.h,
                          volume=8.0,
                          ave_captured_carrier_mass=0.02,
                          ave_static_diele_const=0.03,
                          # down does not exist because of hole capture.
                          e_p_matrix_elements=[e_p_elem])
    assert make_e_p_coupling_h.make() == expected


def test_add_initial_e_p_coupling(e_p_coupling):
    e_p_coupling.e_p_matrix_elements[0].inner_products = []
    spin, kpoint = 0, 1
    wswqs = {(spin, kpoint): {(1, 2): 3.0 + 4.0j}}

    add_inner_products(e_p_coupling, wswqs, dQ=3.0)
    actual = e_p_coupling.e_p_matrix_elements[0].inner_products
    expected = [InnerProduct(inner_product=5.0, dQ=3.0, used_for_fitting=True)]

    assert actual == expected
