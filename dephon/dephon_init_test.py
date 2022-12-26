# -*- coding: utf-8 -*-
#  Copyright (c) 2022 Kumagai group.
from math import sqrt

import pytest
from pymatgen.core import Lattice, Structure, Element
from vise.tests.helpers.assertion import assert_json_roundtrip

from dephon.dephon_init import get_dR, MinimumPointInfo


def test_get_dR():
    lattice = Lattice.orthorhombic(10, 20, 30)
    structure1 = Structure(lattice, ["H"], [[0.0, 0.0, 0.0]])
    structure2 = Structure(lattice, ["H"], [[0.1, 0.1, 0.1]])
    assert get_dR(structure1, structure2) == sqrt(1**2+2**2+3**2)


@pytest.fixture
def minimum_point_info(ground_structure):
    return MinimumPointInfo(charge=-1,
                            structure=ground_structure,
                            bare_energy=-100.0,
                            energy=10.0,
                            correction_energy=1.0,
                            initial_site_symmetry="4mm",
                            final_site_symmetry="2/m",
                            parsed_dir="/path/to/min_point")


def test_json_roundtrip(dephon_init, tmpdir):
    assert_json_roundtrip(dephon_init, tmpdir)


def test_minimum_point_info(minimum_point_info):
    assert minimum_point_info.degeneracy_by_symmetry_reduction == 2


def test_dephon_init_dQ(dephon_init):
    expected = sqrt((0.1*10)**2*6 * Element.H.atomic_mass)
    assert dephon_init.dQ == pytest.approx(expected)


def test_dephon_init_min_point_info_from_charge(dephon_init):
    actual = dephon_init.min_point_info_from_charge(charge=1)
    assert actual == dephon_init.states[1]


def test_dephon_init_volume(dephon_init):
    assert dephon_init.volume == 1000.0


def test_ccd_string(dephon_init):
    actual = dephon_init.__str__()
    print(dephon_init)
    expected = """name: Va_O
vbm             1.000  supercell vbm  1.100
cbm             3.000  supercell cbm  2.900
dQ (amu^0.5 Å)  2.459
dR (Å)          2.449
M (amu)         1.008
------------------------------------------------------------
  q   initial symm     final symm    energy    correction    corrected energy     ZPL
  0       2mm                   2    11.000        -1.000              10.000
  1       2mm                   2    12.000        -1.000              11.000  -1.000"""
    assert actual == expected
