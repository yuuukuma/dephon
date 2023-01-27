# -*- coding: utf-8 -*-
#  Copyright (c) 2022 Kumagai group.
import pytest

from dephon.capture_rate import CaptureRate, make_capture_rate


@pytest.fixture
def capture_rate():
    return CaptureRate(Wif=1.23456,
                       phonon_overlaps=[0.1234, 1.1234, 2.1234],
                       temperatures=[100, 200, 300],
                       degeneracy=4.0,
                       volume=1000.1234)


def test_capture_rate_str(capture_rate):
    print(capture_rate)


def test_make_capture_rate(dephon_init, ccd, e_p_coupling):
    cap_rate = make_capture_rate(dephon_init, ccd, e_p_coupling, [300])
    print(cap_rate)