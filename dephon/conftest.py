# -*- coding: utf-8 -*-
#  Copyright (c) 2022 Kumagai group.
from pathlib import Path

import pytest
from pydefect.analyzer.defect_energy import DefectEnergy
from pymatgen.core import Structure

from dephon.config_coord import ImageStructure, CcdInit


@pytest.fixture(scope="session")
def test_files():
    return Path(__file__).parent / "test_files"


@pytest.fixture(scope="session")
def initial_structure():
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
def final_structure():
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
def ccd_init(initial_structure, final_structure, intermediate_structure):
    return CcdInit(
        name="Va_O",
        initial_structure=initial_structure,
        final_structure=final_structure,
        initial_charge=0,
        final_charge=-1,
        initial_energy=DefectEnergy(-1.0, {"test": -1.0}),
        final_energy=DefectEnergy(-2.0, {"test": -2.0}),
        i_to_f_image_structures=[ImageStructure(intermediate_structure, 0.5)],
        f_to_i_image_structures=[ImageStructure(intermediate_structure, 0.5)])