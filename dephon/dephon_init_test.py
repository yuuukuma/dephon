# -*- coding: utf-8 -*-
#  Copyright (c) 2022 Kumagai group.
from math import sqrt

import pytest
from pymatgen.core import Lattice, Structure, Element
from vise.tests.helpers.assertion import assert_json_roundtrip, \
    assert_dataclass_almost_equal

from dephon.config_coord_test import band_edges
from dephon.dephon_init import get_dR, MinimumPointInfo, DephonInit
from dephon.enum import Carrier


def test_get_dR():
    lattice = Lattice.orthorhombic(10, 20, 30)
    structure1 = Structure(lattice, ["H"], [[0.0, 0.0, 0.0]])
    structure2 = Structure(lattice, ["H"], [[0.1, 0.1, 0.1]])
    assert get_dR(structure1, structure2) == sqrt(1**2+2**2+3**2)


@pytest.fixture
def minimum_point_info(ground_structure):
    return MinimumPointInfo(charge=-1,
                            structure=ground_structure,
                            energy=10.0,
                            energy_correction=1.0,
                            initial_site_symmetry="4mm",
                            final_site_symmetry="2/m",
                            carriers=[Carrier.hole],
                            parsed_dir="/path/to/min_point")


@pytest.fixture
def dephon_init_p_capture(ground_structure, excited_structure):
    ground_state = MinimumPointInfo(charge=0,
                                    structure=ground_structure,
                                    energy=11.0,
                                    energy_correction=-1.0,
                                    initial_site_symmetry="2mm",
                                    final_site_symmetry="2",
                                    carriers=[],
                                    parsed_dir="/path/to/ground")
    excited_state = MinimumPointInfo(charge=-1,
                                     structure=excited_structure,
                                     energy=12.0,
                                     energy_correction=-1.0,
                                     initial_site_symmetry="2mm",
                                     final_site_symmetry="2",
                                     carriers=[Carrier.hole],
                                     parsed_dir="/path/to/excited")
    # transition level = 1.0 from VBM
    return DephonInit(name="Va_O", ground_state=ground_state,
                      excited_state=excited_state, **band_edges)


@pytest.fixture
def dephon_init_n_capture(ground_structure, excited_structure):
    ground_state = MinimumPointInfo(charge=0,
                                    structure=ground_structure,
                                    energy=11.0,
                                    energy_correction=-1.0,
                                    initial_site_symmetry="2mm",
                                    final_site_symmetry="2",
                                    carriers=[],
                                    parsed_dir="/path/to/ground")
    excited_state = MinimumPointInfo(charge=1,
                                     structure=excited_structure,
                                     energy=12.0,
                                     energy_correction=-1.0,
                                     initial_site_symmetry="2mm",
                                     final_site_symmetry="2",
                                     carriers=[Carrier.electron],
                                     parsed_dir="/path/to/excited")
    # transition level = -1.0 from CBM
    return DephonInit(name="Va_O", ground_state=ground_state,
                      excited_state=excited_state, **band_edges)


def test_json_roundtrip(dephon_init_p_capture, tmpdir):
    assert_json_roundtrip(dephon_init_p_capture, tmpdir)


def test_minimum_point_info(minimum_point_info):
    assert minimum_point_info.degeneracy_by_symmetry_reduction == 2


def test_dephon_init_ground_state_w_pn(dephon_init_p_capture, ground_structure):
    actual = dephon_init_p_capture.ground_state_w_pn
    expected = MinimumPointInfo(charge=0,
                                structure=ground_structure,
                                energy=12.0,
                                energy_correction=-1.0,
                                initial_site_symmetry="2mm",
                                final_site_symmetry="2",
                                carriers=[Carrier.hole, Carrier.electron],
                                parsed_dir="/path/to/ground")
    assert_dataclass_almost_equal(actual, expected)


def test_dephon_init_dQ(dephon_init_p_capture):
    expected = sqrt((0.1*10)**2*6 * Element.H.atomic_mass)
    assert dephon_init_p_capture.dQ == pytest.approx(expected)


def test_dephon_init_semiconductor_type(dephon_init_p_capture,
                                        dephon_init_n_capture):
    assert dephon_init_p_capture.semiconductor_type == "p"
    assert dephon_init_n_capture.semiconductor_type == "n"


def test_dephon_init_volume(dephon_init_p_capture):
    assert dephon_init_p_capture.volume == 1000.0


def test_ccd_string(dephon_init_p_capture):
    actual = dephon_init_p_capture.__str__()
    print(dephon_init_p_capture)
    expected = """name: Va_O
semiconductor type:  p-type
vbm             1.000  supercell vbm  1.100
cbm             2.000  supercell cbm  1.900
dQ (amu^0.5 Å)  2.459
dR (Å)          2.449
M (amu)         1.008
------------------------------------------------------------
  q   carrier    initial symm     final symm    energy    correction    corrected energy    ZPL
  0     h+e          2mm                   2    12.000        -1.000              11.000
 -1      h           2mm                   2    12.000        -1.000              11.000  0.000
  0                  2mm                   2    11.000        -1.000              10.000  1.000"""
    assert actual == expected
