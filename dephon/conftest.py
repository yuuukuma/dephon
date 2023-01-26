# -*- coding: utf-8 -*-
#  Copyright (c) 2022 Kumagai group.
from copy import copy
from pathlib import Path

import pytest
from pydefect.analyzer.band_edge_states import LocalizedOrbital
from pymatgen.core import Structure, Lattice
from pymatgen.electronic_structure.core import Spin

from dephon.config_coord import SinglePointInfo, Ccd, SingleCcd, SingleCcdId
from dephon.dephon_init import MinimumPointInfo, DephonInit, NearEdgeState
from dephon.ele_phon_coupling import InnerProduct, EPMatrixElement, EPCoupling
from dephon.enum import Carrier


@pytest.fixture(scope="session")
def test_files():
    return Path(__file__).parent / "test_files"


@pytest.fixture(scope="session")
def excited_structure():
    return Structure.from_str("""H He
10.0
1.0 0.0 0.0
0.0 1.0 0.0
0.0 0.0 1.0
H
6
direct 
0.1 0.0 0.0
0.9 0.0 0.0
0.0 0.1 0.0
0.0 0.9 0.0
0.0 0.0 0.1
0.0 0.0 0.9""", fmt="poscar")


@pytest.fixture(scope="session")
def ground_structure():
    return Structure.from_str("""Mg4 O4
10.0
1.0 0.0 0.0
0.0 1.0 0.0
0.0 0.0 1.0
H
6
direct
0.2 0.0 0.0
0.8 0.0 0.0
0.0 0.2 0.0
0.0 0.8 0.0
0.0 0.0 0.2
0.0 0.0 0.8""", fmt="poscar")


@pytest.fixture(scope="session")
def intermediate_structure():
    return Structure.from_str("""Mg4 O4
10.0
1.0 0.0 0.0
0.0 1.0 0.0
0.0 0.0 1.0
H
6
direct
0.15 0.0 0.0
0.85 0.0 0.0
0.0 0.15 0.0
0.0 0.85 0.0
0.0 0.0 0.15
0.0 0.0 0.85""", fmt="poscar")


@pytest.fixture
def dephon_init(ground_structure, excited_structure):
    orb_info = LocalizedOrbital(band_idx=2,
                                ave_energy=2.0,
                                occupation=1.0,
                                orbitals={"O": [0.0, 1.0, 0.0]})
    vb = NearEdgeState(band_index=1,
                       kpt_coord=[0.0]*3,
                       kpt_index=1,
                       kpt_weight=1.0,
                       eigenvalue=1.0,
                       occupation=1.0)
    cb = NearEdgeState(band_index=2,
                       kpt_coord=[0.0]*3,
                       kpt_index=1,
                       kpt_weight=1.0,
                       eigenvalue=3.0,
                       occupation=0.0)
    cb_w_lo = copy(cb)
    cb_w_lo.band_index = 3

    va_o1_0 = MinimumPointInfo(charge=0,
                               structure=ground_structure,
                               energy=11.0,
                               correction_energy=-1.0,
                               magnetization=0.0,
                               localized_orbitals=[[], []],
                               initial_site_symmetry="2mm",
                               final_site_symmetry="2",
                               parsed_dir="/path/to/Va_O1_0",
                               valence_bands=[[vb], [vb]],
                               conduction_bands=[[cb], [cb]])
    va_o1_1 = MinimumPointInfo(charge=1,
                               structure=excited_structure,
                               energy=12.0,
                               correction_energy=-1.0,
                               magnetization=1.0,
                               localized_orbitals=[[], [orb_info]],
                               initial_site_symmetry="2mm",
                               final_site_symmetry="2",
                               parsed_dir="/path/to/Va_O1_1",
                               valence_bands=[[vb], [vb]],
                               conduction_bands=[[cb], [cb_w_lo]])
    # transition level = -1.0 from CBM
    return DephonInit(defect_name="Va_O", states=[va_o1_0, va_o1_1],
                      vbm=1.0, cbm=3.0, supercell_vbm=1.1, supercell_cbm=2.9,
                      ave_electron_mass=11.0, ave_hole_mass=12.0,
                      ave_static_diele_const=13.0)


@pytest.fixture
def ccd(excited_structure, ground_structure, intermediate_structure):
    return Ccd(defect_name="test",
               ccds=[
                   SingleCcd(SingleCcdId(name="excited"), charge=0,
                             point_infos=[SinglePointInfo(-1.0, -0.1, 2.1, False, used_for_fitting=True),
                                          SinglePointInfo(0.0, 0.0, 1.1, False, used_for_fitting=True),
                                          SinglePointInfo(1.0, 0.1, 2.2, False, used_for_fitting=True)]),
                   SingleCcd(SingleCcdId(name="ground"), charge=1,
                             point_infos=[SinglePointInfo(-1.0, -0.1, 1.1, False, used_for_fitting=False),
                                          SinglePointInfo(0.0, 0.0, 0.1, False, used_for_fitting=True),
                                          SinglePointInfo(1.0, 0.1, 1.2, False, used_for_fitting=True)])])


@pytest.fixture(scope="session")
def sc_structure():
    lattice = Lattice.cubic(1.0)
    coords = [[0.0, 0.0, 0.0]]
    return Structure(lattice=lattice, species=["H"], coords=coords)


@pytest.fixture
def e_p_matrix_elem():
    ip_1 = InnerProduct(inner_product=20.0, dQ=-1.0, used_for_fitting=False)
    ip_2 = InnerProduct(inner_product=1.0, dQ=0.0, used_for_fitting=True)
    ip_3 = InnerProduct(inner_product=2.0, dQ=1.0, used_for_fitting=True)

    return EPMatrixElement(band_edge_index=1,
                           defect_band_index=2,
                           spin=Spin.down,
                           eigenvalue_diff=0.1,
                           kpt_idx=1,
                           kpt_coord=[0.0, 0.0, 0.0],
                           inner_products=[ip_1, ip_2, ip_3])


@pytest.fixture
def e_p_coupling(e_p_matrix_elem):
    return EPCoupling(
        charge=1,
        disp=0.0,
        captured_carrier=Carrier.e,
        volume=100.0,
        ave_captured_carrier_mass=1.0,
        ave_static_diele_const=2.0,
        e_p_matrix_elements=[e_p_matrix_elem])
