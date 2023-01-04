# -*- coding: utf-8 -*-
#  Copyright (c) 2022 Kumagai group.
from pathlib import Path

import pytest
from pymatgen.core import Structure, Lattice

from dephon.config_coord import SinglePointInfo, Ccd, SingleCcd
from dephon.dephon_init import MinimumPointInfo, DephonInit


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


band_edges = dict(vbm=1.0, cbm=3.0, supercell_vbm=1.1, supercell_cbm=2.9)


@pytest.fixture
def dephon_init(ground_structure, excited_structure):
    va_o1_0 = MinimumPointInfo(charge=0,
                               structure=ground_structure,
                               energy=11.0,
                               correction_energy=-1.0,
                               initial_site_symmetry="2mm",
                               final_site_symmetry="2",
                               parsed_dir="/path/to/Va_O1_0")
    va_o1_1 = MinimumPointInfo(charge=1,
                               structure=excited_structure,
                               energy=12.0,
                               correction_energy=-1.0,
                               initial_site_symmetry="2mm",
                               final_site_symmetry="2",
                               parsed_dir="/path/to/Va_O1_1")
    # transition level = -1.0 from CBM
    return DephonInit(defect_name="Va_O", states=[va_o1_0, va_o1_1],
                      **band_edges)


@pytest.fixture
def ccd(excited_structure, ground_structure, intermediate_structure):
    return Ccd(defect_name="test",
               ccds=[
                   SingleCcd(name="excited",
                             charge=0,
                             point_infos=[SinglePointInfo(-1.0, -0.1, 2.1, False, used_for_fitting=True),
                                          SinglePointInfo(0.0, 0.0, 1.1, False, used_for_fitting=True),
                                          SinglePointInfo(1.0, 0.1, 2.2, False, used_for_fitting=True)]),
                   SingleCcd(name="ground",
                             charge=1,
                             point_infos=[SinglePointInfo(-1.0, -0.1, 1.1, False, used_for_fitting=False),
                                          SinglePointInfo(0.0, 0.0, 0.1, False, used_for_fitting=True),
                                          SinglePointInfo(1.0, 0.1, 1.2, False, used_for_fitting=True)])])


@pytest.fixture(scope="session")
def sc_structure():
    lattice = Lattice.cubic(1.0)
    coords = [[0.0, 0.0, 0.0]]
    return Structure(lattice=lattice, species=["H"], coords=coords)