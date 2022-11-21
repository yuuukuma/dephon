# -*- coding: utf-8 -*-
#  Copyright (c) 2022 Kumagai group.
from math import sqrt

import numpy as np
import pytest
from pymatgen.core import Element
from vise.tests.helpers.assertion import assert_json_roundtrip

from dephon.config_coord import Ccd, ImageStructureInfo, ccd_plt


def test_ccd_init_to_json_roundtrip(ccd_init, tmpdir):
    assert_json_roundtrip(ccd_init, tmpdir)


def test_ccd_to_json_roundtrip(ccd, tmpdir):
    assert_json_roundtrip(ccd, tmpdir)


def test_ccd_dQ(ccd_init):
    expected = sqrt((0.1*10)**2*6 * Element.H.atomic_mass)
    assert ccd_init.dQ == pytest.approx(expected)


def test_ccd_string(ccd_init):
    actual = ccd_init.__str__()
    expected = """Excited state:  Va_O_0 + e-  energy:  -1  correction:  -1
Ground state:   Va_O_-1      energy:  -2  correction:  -2"""
    assert actual == expected


@pytest.fixture
def ccd():
    return Ccd(dQ=10.0,
               excited_image_infos=[
                   ImageStructureInfo(-0.2, 11.4),
                   ImageStructureInfo(0.0, 11.0),
                   ImageStructureInfo(0.2, 11.4),
                   ImageStructureInfo(0.4, 12.6),
                   ImageStructureInfo(0.6, 14.6),
                   ImageStructureInfo(0.8, 17.4),
                   ImageStructureInfo(1.0, 21.0)],
               ground_image_infos=[
                   ImageStructureInfo(-0.2, 10.4),
                   ImageStructureInfo(0.0, 10.0),
                   ImageStructureInfo(0.2, 10.4),
                   ImageStructureInfo(0.4, 11.6),
                   ImageStructureInfo(0.6, 13.6),
                   ImageStructureInfo(0.8, 16.4),
                   ImageStructureInfo(1.0, 20.0)])


def test_ccd(ccd):
    assert ccd.ground_dQs == [-2.0, 0.0, 2.0, 4.0, 6.0, 8.0, 10.0]
    np.testing.assert_array_almost_equal(ccd.excited_dQs, [12.0, 10.0, 8.0, 6.0, 4.0, 2.0, 0.0])


def test_plot_ccd(ccd):
    plt = ccd_plt(ccd)
    plt.show()

""" 
TODO: 
1. add defect position
2. consider how to handle the small difference of origin.
"""