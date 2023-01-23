# -*- coding: utf-8 -*-
#  Copyright (c) 2022 Kumagai group.
from copy import copy
from math import sqrt

import pytest
from pydefect.analyzer.band_edge_states import LocalizedOrbital
from pymatgen.core import Lattice, Structure, Element
from pymatgen.electronic_structure.core import Spin
from vise.tests.helpers.assertion import assert_json_roundtrip

from dephon.dephon_init import get_dR, MinimumPointInfo, NearEdgeState
from dephon.enum import Carrier


def test_get_dR():
    lattice = Lattice.orthorhombic(10, 20, 30)
    structure1 = Structure(lattice, ["H"], [[0.0, 0.0, 0.0]])
    structure2 = Structure(lattice, ["H"], [[0.1, 0.1, 0.1]])
    assert get_dR(structure1, structure2) == sqrt(1**2+2**2+3**2)


vb_nes = NearEdgeState(band_index=1,
                       kpt_coord=[0.0]*3,
                       kpt_weight=1.0,
                       kpt_index=1,
                       eigenvalue=1.5,
                       occupation=1.0)


@pytest.fixture
def minimum_point_info(ground_structure):
    orb_info = LocalizedOrbital(band_idx=2,
                                ave_energy=2.0,
                                occupation=1.0,
                                orbitals={"O": [0.0, 1.0, 0.0]})
    cb_nes_up = NearEdgeState(band_index=2,
                              kpt_coord=[0.0]*3,
                              kpt_weight=1.0,
                              kpt_index=1,
                              eigenvalue=2.5,
                              occupation=0.0)
    cb_nes_down = copy(cb_nes_up)
    cb_nes_down.band_index = 3
    return MinimumPointInfo(charge=-1,
                            structure=ground_structure,
                            energy=-90.0,
                            correction_energy=1.0,
                            magnetization=1.0,
                            localized_orbitals=[[], [orb_info]],
                            valence_bands=[[vb_nes]],
                            conduction_bands=[[cb_nes_up], [cb_nes_down]],
                            initial_site_symmetry="4mm",
                            final_site_symmetry="2/m",
                            parsed_dir="/path/to/min_point")


def test_json_roundtrip(dephon_init, tmpdir):
    assert_json_roundtrip(dephon_init, tmpdir)


def test_minimum_point_info_degeneracy_by_symm_reduction(minimum_point_info):
    assert minimum_point_info.degeneracy_by_symmetry_reduction == 2


def test_minimum_point_info_near_edge_states(minimum_point_info):
    actual = minimum_point_info.near_edge_states(
        captured_carrier=Carrier.h, spin=Spin.down)
    assert actual == [vb_nes]


def test_minimum_point_info_relevant_band_indices(minimum_point_info):
    expected = {(1, 1): [1, 2], (2, 1): [1, 3, 2]}
    assert minimum_point_info.relevant_band_indices == expected


def test_dephon_init_dQ(dephon_init):
    expected = sqrt((0.1*10)**2*6 * Element.H.atomic_mass)
    assert dephon_init.dQ == pytest.approx(expected)


def test_dephon_init_min_point_info_from_charge(dephon_init):
    actual = dephon_init.min_info_from_charge(charge=1)
    assert actual == dephon_init.states[1]


def test_dephon_init_volume(dephon_init):
    assert dephon_init.volume == 1000.0


def test_ccd_string(dephon_init):
    actual = dephon_init.__str__()
    print(dephon_init)
    expected = """name: Va_O
vbm                  1.000  supercell vbm  1.100
cbm                  3.000  supercell cbm  2.900
dQ (amu^0.5 Å)       2.459
dR (Å)               2.449
M (amu)              1.008
electron mass (m0)  11.000
hole mass (m0)      12.000
static diele        13.000
------------------------------------------------------------
  q   ini symm     final symm    energy    correction    corrected energy    magnetization   localized state idx      ZPL
  0     2mm                 2    11.000        -1.000              10.000            0.000
  1     2mm                 2    12.000        -1.000              11.000            1.000         down-2          -1.000"""
    assert actual == expected
