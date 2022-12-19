# -*- coding: utf-8 -*-
#  Copyright (c) 2022 Kumagai group.
from math import sqrt

import pytest
from pymatgen.core import Element, Structure, Lattice
from vise.tests.helpers.assertion import assert_json_roundtrip, \
    assert_dataclass_almost_equal

from dephon.config_coord import Ccd, ImageStructureInfo, CcdPlotter, \
    get_dR, CcdInit, MinimumPointInfo
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


band_edges = dict(vbm=1.0, cbm=2.0, supercell_vbm=1.1, supercell_cbm=1.9)


@pytest.fixture
def ccd_init_p_capture(ground_structure, excited_structure):
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
    return CcdInit(name="Va_O", ground_state=ground_state,
                   excited_state=excited_state, **band_edges)


@pytest.fixture
def ccd_init_n_capture(ground_structure, excited_structure):
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
    return CcdInit(name="Va_O", ground_state=ground_state,
                   excited_state=excited_state, **band_edges)


def test_json_roundtrip(ccd_init_p_capture, tmpdir):
    assert_json_roundtrip(ccd_init_p_capture, tmpdir)


def test_minimum_point_info(minimum_point_info):
    assert minimum_point_info.degeneracy_by_symmetry_reduction == 2


def test_ccd_init_ground_state_w_pn(ccd_init_p_capture, ground_structure):
    actual = ccd_init_p_capture.ground_state_w_pn
    expected = MinimumPointInfo(charge=0,
                                structure=ground_structure,
                                energy=12.0,
                                energy_correction=-1.0,
                                initial_site_symmetry="2mm",
                                final_site_symmetry="2",
                                carriers=[Carrier.hole, Carrier.electron],
                                parsed_dir="/path/to/ground")
    assert_dataclass_almost_equal(actual, expected)


def test_ccd_init_dQ(ccd_init_p_capture):
    expected = sqrt((0.1*10)**2*6 * Element.H.atomic_mass)
    assert ccd_init_p_capture.dQ == pytest.approx(expected)


def test_ccd_init_minority_carrier(ccd_init_p_capture, ccd_init_n_capture):
    assert ccd_init_p_capture.minority_carrier == Carrier.hole
    assert ccd_init_n_capture.minority_carrier == Carrier.electron


def test_ccd_init_semiconductor_type(ccd_init_p_capture, ccd_init_n_capture):
    assert ccd_init_p_capture.semiconductor_type == "p-type"
    assert ccd_init_n_capture.semiconductor_type == "n-type"


def test_ccd_string(ccd_init_p_capture):
    actual = ccd_init_p_capture.__str__()
    print(ccd_init_p_capture)
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


def test_image_structure_info():
    imag_structure_info = ImageStructureInfo(dQ=1.0, energy=2.0, correction=3.0)
    assert imag_structure_info.corrected_energy == 2.0 + 3.0


@pytest.fixture
def ccd():
    return Ccd(name="Va_O1",
               image_infos={"q=0": [ImageStructureInfo(9.4, 1.0, 12.0),
                                    ImageStructureInfo(10.0, 1.0, 10.0),
                                    ImageStructureInfo(10.4, 1.0, 6.0),
                                    ImageStructureInfo(10.8, 1.0, 4.0),
                                    ImageStructureInfo(13.6, 1.0, 3.0),
                                    ImageStructureInfo(16.4, 1.0, 2.0),
                                    ImageStructureInfo(20.0, 1.0, 0.0)],
                            "q=1": [ImageStructureInfo(8.4, 1.0, -2.0),
                                    ImageStructureInfo(9.0, 1.0, -3.5),
                                    ImageStructureInfo(9.4, 1.0, -2.0),
                                    ImageStructureInfo(10.6, 1.0, -1.0),
                                    ImageStructureInfo(12.6, 1.0, 1.0),
                                    ImageStructureInfo(15.4, 1.0, 4.0),
                                    ImageStructureInfo(19.0, 1.0, 10.0)]},
               fitting_q_ranges={"q=1": [10.0, 20.0]})


def test_ccd_omega_(ccd):
    ccd.set_q_range("q=0", 9.0, 21.0)
    assert ccd.omega(image_name="q=0") == 0.033121755227466375


def test_ccd_omega_str(ccd):
    print(ccd)


def test_ccd_lowest_energy(ccd):
    assert ccd.lowest_energy == -2.5


def test_plot_ccd(ccd):
    plotter = CcdPlotter(ccd)
    plotter.construct_plot()
    plotter.plt.show()

""" 
TODO: 
plot

1. add defect position
2. consider how to handle the small difference of origin.
"""