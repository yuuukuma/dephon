# -*- coding: utf-8 -*-
#  Copyright (c) 2023 Kumagai group.
import numpy as np
import pytest
from vise.tests.helpers.assertion import assert_json_roundtrip

from dephon.ccd.screening_charge import ScreeningCharge


@pytest.fixture
def screening_charge():
    return ScreeningCharge(1, -1, np.eye(3), np.diag([2, 3, 4]), 0.1)


def test_screening_charge(screening_charge, tmpdir):
    assert_json_roundtrip(screening_charge, tmpdir)


def test_properties(screening_charge):
    e_prime = 1.0 - 1 / 4
    assert screening_charge.charge == -1 * e_prime + 1 * e_prime * 0.1