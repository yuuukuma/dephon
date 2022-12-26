# -*- coding: utf-8 -*-
#  Copyright (c) 2022 Kumagai group.
import pytest
from vise.tests.helpers.assertion import assert_msonable, assert_yaml_roundtrip

from dephon.corrections import DephonCorrection
from dephon.enum import CorrectionType


@pytest.fixture
def dephon_correction():
    return DephonCorrection(energy=1.0, correction_type=CorrectionType.extended_FNV)


def test_correction_msonable(dephon_correction):
    assert_msonable(dephon_correction)


def test_correction_yaml_round_trip(dephon_correction, tmpdir):
    expected_str = f"""correction_type: extended FNV
energy: 1.0
"""
    assert_yaml_roundtrip(dephon_correction, tmpdir, expected_str)

