# -*- coding: utf-8 -*-
#  Copyright (c) 2022 Kumagai group.
from pathlib import Path

import pytest
from pymatgen.core import Structure, Lattice

from dephon.config_coord import SinglePointInfo, Ccd, SingleCcd
from dephon.enum import CorrectionEnergyType


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
def ccd(excited_structure, ground_structure, intermediate_structure):
    return Ccd(defect_name="test",
               correction_energy_type=CorrectionEnergyType.extended_FNV,
               image_infos_list=[
                   SingleCcd("excited",
                             [SinglePointInfo(0.0, 0.1, 0.1, False, -3.0),
                              SinglePointInfo(1.0, 0.2, 1.1, False, -3.0)]),
                   SingleCcd("ground",
                             [SinglePointInfo(0.0, 0.9, 0.2, False, -4.0),
                              SinglePointInfo(1.0, 0.8, 1.2, False, -4.0)])])


@pytest.fixture(scope="session")
def sc_structure():
    lattice = Lattice.cubic(1.0)
    coords = [[0.0, 0.0, 0.0]]
    return Structure(lattice=lattice, species=["H"], coords=coords)