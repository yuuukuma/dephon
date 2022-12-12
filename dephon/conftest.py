# -*- coding: utf-8 -*-
#  Copyright (c) 2022 Kumagai group.
from pathlib import Path

import pytest
from pydefect.analyzer.defect_energy import DefectEnergy
from pymatgen.core import Structure, Lattice

from dephon.config_coord import ImageStructureInfo, CcdInit, Ccd


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
def ccd_init(excited_structure, ground_structure, intermediate_structure):
    return CcdInit(
        name="Va_O",
        excited_structure=excited_structure,
        ground_structure=ground_structure,
        excited_charge=0,
        ground_charge=-1,
        excited_energy=-1.0,
        excited_energy_correction=-1.0,
        ground_energy=-2.0,
        ground_energy_correction=-2.0,
        vbm=1.0,
        cbm=2.0,
        supercell_vbm=1.1,
        supercell_cbm=1.9)


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