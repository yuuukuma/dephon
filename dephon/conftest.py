# -*- coding: utf-8 -*-
#  Copyright (c) 2022 Kumagai group.
from pathlib import Path

import pytest
from pymatgen.core import Structure, Lattice

from dephon.config_coord import ImageStructureInfo, CcdInit, Ccd, \
    MinimumPointInfo


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
def minimum_point_info(ground_structure):
    return MinimumPointInfo(charge=-1,
                            structure=ground_structure,
                            energy=-2.0123456,
                            energy_correction=-2.1234567,
                            initial_site_symm="4mm",
                            final_site_symm="2/m")


@pytest.fixture
def ccd_init(minimum_point_info, excited_structure):
    excited_state = MinimumPointInfo(charge=0, structure=excited_structure,
                                     energy=-1.0123456,
                                     energy_correction=-1.1234567,
                                     initial_site_symm="2mm",
                                     final_site_symm="2")

    return CcdInit(
        name="Va_O",
        ground_state=minimum_point_info,
        excited_state=excited_state,
        vbm=1.0123456,
        cbm=2.0123456,
        supercell_vbm=1.1023456,
        supercell_cbm=1.9123456)


@pytest.fixture
def ccd(excited_structure, ground_structure, intermediate_structure):
    return Ccd(image_infos={
        "excited": [ImageStructureInfo(0.0, 0.1, -3.0),
                    ImageStructureInfo(1.0, 1.1, -3.0)],
        "ground": [ImageStructureInfo(0.0, 0.2, -4.0),
                   ImageStructureInfo(1.0, 1.2, -4.0)]},
        name="test")


@pytest.fixture(scope="session")
def sc_structure():
    lattice = Lattice.cubic(1.0)
    coords = [[0.0, 0.0, 0.0]]
    return Structure(lattice=lattice, species=["H"], coords=coords)