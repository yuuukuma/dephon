# -*- coding: utf-8 -*-
#  Copyright (c) 2022 Kumagai group.
from math import sqrt

import pytest
from pymatgen.core import Element, Structure, Lattice
from vise.tests.helpers.assertion import assert_json_roundtrip

from dephon.config_coord import Ccd, ImageStructureInfo, CcdPlotter, \
    get_dR


def test_get_dR():
    lattice = Lattice.orthorhombic(10, 20, 30)
    structure1 = Structure(lattice, ["H"], [[0.0, 0.0, 0.0]])
    structure2 = Structure(lattice, ["H"], [[0.1, 0.1, 0.1]])
    assert get_dR(structure1, structure2) == sqrt(1**2+2**2+3**2)


def test_ccd_init_to_json_roundtrip(ccd_init, tmpdir):
    assert_json_roundtrip(ccd_init, tmpdir)


def test_ccd_to_json_roundtrip(ccd, tmpdir):
    assert_json_roundtrip(ccd, tmpdir)


def test_ccd_dQ(ccd_init):
    expected = sqrt((0.1*10)**2*6 * Element.H.atomic_mass)
    assert ccd_init.dQ == pytest.approx(expected)


def test_ccd_string(ccd_init):
    actual = ccd_init.__str__()
    print(ccd_init)
    expected = """name: Va_O
vbm             1.012  supercell vbm  1.102
cbm             2.012  supercell cbm  1.912
dQ (amu^0.5 Å)  2.459
dR (Å)          2.449
M (amu)         1.008
------------------------------------------------------------
 state    initial symm    final symm     energy    correction    corrected energy
Va_O_0        2mm             2          -1.012        -1.123              -2.136
Va_O_-1       2mm             2m         -2.012        -2.123              -4.136
ZPL: 2.000"""
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