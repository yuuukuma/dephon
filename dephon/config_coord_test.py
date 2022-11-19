# -*- coding: utf-8 -*-
#  Copyright (c) 2022 Kumagai group.
from math import sqrt

import pytest
from pymatgen.core import Element
from vise.tests.helpers.assertion import assert_json_roundtrip


def test_ccd_to_json_roundtrip(ccd_init, tmpdir):
    assert_json_roundtrip(ccd_init, tmpdir)


def test_ccd_dQ(ccd_init):
    expected = sqrt((0.1*10)**2*6 * Element.H.atomic_mass)
    assert ccd_init.dQ == pytest.approx(expected)


def test_ccd_string(ccd_init):
    actual = ccd_init.__str__()
    expected = """Name: Va_O
transition: Va_O_0 + e- -> Va_O_-1"""
    assert actual == expected


""" 
TODO: 
1. add defect position
2. consider how to handle the small difference of origion.
"""