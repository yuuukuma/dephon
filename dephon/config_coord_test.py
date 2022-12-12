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
    expected = """dQ              2.46
dR              2.45
M               1.01
Excited state:  Va_O_0 + e-  energy:  -1  correction:  -1
Ground state:   Va_O_-1      energy:  -2  correction:  -2"""
    assert actual == expected


def test_image_structure_info():
    imag_structure_info = ImageStructureInfo(dQ=1.0, energy=2.0, correction=3.0)
    assert imag_structure_info.corrected_energy == 2.0 + 3.0


@pytest.fixture
def ccd():
    return Ccd(name="Va_O1",
               image_infos={"q=0": [ImageStructureInfo(9.4, 1.0, 12.0),
                                    ImageStructureInfo(10.0, 1.0, 10.0),
                                    ImageStructureInfo(10.4, 1.0, 8.0),
                                    ImageStructureInfo(11.6, 1.0, 6.0),
                                    ImageStructureInfo(13.6, 1.0, 4.0),
                                    ImageStructureInfo(16.4, 1.0, 2.0),
                                    ImageStructureInfo(20.0, 1.0, 0.0)],
                            "q=1": [ImageStructureInfo(8.4, 1.0, -2.0),
                                    ImageStructureInfo(9.0, 1.0, 0.0),
                                    ImageStructureInfo(9.4, 1.0, 2.0),
                                    ImageStructureInfo(10.6, 1.0, 4.0),
                                    ImageStructureInfo(12.6, 1.0, 6.0),
                                    ImageStructureInfo(15.4, 1.0, 8.0),
                                    ImageStructureInfo(19.0, 1.0, 10.0)]})


def test_ccd_min_energy(ccd):
    assert ccd.min_energy == -1.0


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