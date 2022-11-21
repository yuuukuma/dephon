# -*- coding: utf-8 -*-
#  Copyright (c) 2022 Kumagai group.
from math import sqrt

import pytest
from pymatgen.core import Element
from vise.tests.helpers.assertion import assert_json_roundtrip


def test_ccd_init_to_json_roundtrip(ccd_init, tmpdir):
    assert_json_roundtrip(ccd_init, tmpdir)


def test_ccd_to_json_roundtrip(ccd, tmpdir):
    assert_json_roundtrip(ccd, tmpdir)


def test_ccd_dQ(ccd_init):
    expected = sqrt((0.1*10)**2*6 * Element.H.atomic_mass)
    assert ccd_init.dQ == pytest.approx(expected)


def test_ccd_string(ccd_init):
    actual = ccd_init.__str__()
    expected = """Name:           Va_O
Excited state:  Va_O_0 + e-  energy:  -1  total correction:  -1  is shallow:
Ground state:   Va_O_-1      energy:  -2  total correction:  -2  is shallow:"""
    assert actual == expected


""" 
TODO: 
1. add defect position
2. consider how to handle the small difference of origin.
"""