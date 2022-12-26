# -*- coding: utf-8 -*-
#  Copyright (c) 2022 Kumagai group.
from copy import deepcopy

import numpy as np
import pytest
from vise.tests.helpers.assertion import assert_dataclass_almost_equal

from dephon.config_coord import Ccd, SinglePointInfo, CcdPlotter, \
    SingleCcd, spline3
from dephon.enum import CorrectionType, Carrier


@pytest.fixture
def single_point_min_info():
    return SinglePointInfo(dQ=1.0, disp_ratio=0.1)


@pytest.fixture
def single_point_max_info():
    return SinglePointInfo(dQ=1.0,
                           disp_ratio=0.1,
                           corrected_energy=5.0,
                           is_shallow=True,
                           correction_method=CorrectionType.extended_FNV,
                           used_for_fitting=True,
                           base_energy=1.0)


def test_single_point_info_corrected_energy(single_point_min_info,
                                            single_point_max_info):
    assert single_point_min_info.relative_energy is None
    assert single_point_max_info.relative_energy == 5.0 - 1.0


def test_image_structure_str(single_point_min_info,
                             single_point_max_info):
    actual = single_point_min_info.__str__()
    expected = """  dQ    disp ratio  corr. energy    relative energy    used for fitting?    is shallow?
1.00          0.10  -               -                  -                    -"""
    assert actual == expected

    actual = single_point_max_info.__str__()
    expected = """  dQ    disp ratio    corr. energy    relative energy  used for fitting?    is shallow?
1.00          0.10            5.00               4.00  True                 True"""
    assert actual == expected


@pytest.fixture
def single_ccd():
    return SingleCcd(
        name="from_0_to_1",
        charge=0,
        carriers=[Carrier.hole, Carrier.electron],
        point_infos=[
            SinglePointInfo(2., 1.0, 3.3, False, used_for_fitting=True),
            SinglePointInfo(1., 0.5, 2.2, False, used_for_fitting=True),
            SinglePointInfo(0., 0.0, 3.4, False, used_for_fitting=True),
            SinglePointInfo(3., 1.5, 3.5, False, used_for_fitting=False)])


def test_single_ccd_sort_single_point_infos(single_ccd):
    actual = single_ccd.point_infos[0]
    expected = SinglePointInfo(0., 0.0, 3.4, False, used_for_fitting=True)
    assert actual == expected


def test_single_ccd_dQs_and_energies(single_ccd):
    single_ccd.set_base_energy(0.1)

    actual = single_ccd.dQs_and_energies(False)
    dQs = [0.0, 1.0, 2.0, 3.0]
    energies = [3.3, 2.1, 3.2, 3.4]
    expected = (dQs, energies)
    assert np.array(actual) == pytest.approx(np.array(expected))

    single_ccd.set_base_energy()

    actual = single_ccd.dQs_and_energies(only_used_for_fitting=True)
    dQs = [0.0, 1.0, 2.0]
    energies = [0.0, -1.2, -0.1]
    expected = (dQs, energies)
    assert np.array(actual) == pytest.approx(np.array(expected))


def test_single_ccd_omega(single_ccd):
    assert single_ccd.omega() == pytest.approx(0.09805287531186566)


def test_energy_shifted_single_ccd(single_ccd):
    actual = single_ccd.energy_shifted_single_ccd(dQ_0_energy=1.0)
    expected = deepcopy(single_ccd)
    for i in expected.point_infos:
        i.base_energy = 2.4
    assert_dataclass_almost_equal(actual, expected)


def test_dQ_reverted_single_ccd(single_ccd):
    actual = single_ccd.dQ_reverted_single_ccd()
    expected = deepcopy(single_ccd)
    expected.point_infos = [
        SinglePointInfo(-1.0, 1.5, 3.5, False, used_for_fitting=False),
        SinglePointInfo(0.0, 1.0, 3.3, False, used_for_fitting=True),
        SinglePointInfo(1.0, 0.5, 2.2, False, used_for_fitting=True),
        SinglePointInfo(2.0, 0.0, 3.4, False, used_for_fitting=True)]
    assert_dataclass_almost_equal(actual, expected)


def test_single_ccd_str(single_ccd):
    single_ccd.set_base_energy()
    actual = single_ccd.__str__()
    print(actual)
    expected = """name: from_0_to_1
omega: 0.10
carriers: h e
  dQ    disp ratio    corr. energy    relative energy  used for fitting?    is shallow?
0.00          0.00            3.40               0.00  True                 False
1.00          0.50            2.20              -1.20  True                 False
2.00          1.00            3.30              -0.10  True                 False
3.00          1.50            3.50               0.10  False                False"""
    assert actual == expected


@pytest.fixture
def ccd(single_ccd):
    excited_state = SingleCcd(
        name="from_1_to_0",
        charge=1,
        point_infos=[SinglePointInfo(3., 1.0, 10.1, False),
                     SinglePointInfo(2., 0.9, 10.2, False),
                     SinglePointInfo(1., 0.8, 10.3, False)])
    return Ccd(defect_name="Va_O1", ccds=[single_ccd, excited_state])


def test_ccd_single_ccd(ccd, single_ccd):
    actual = ccd.single_ccd("from_0_to_1")
    expected = single_ccd
    assert actual == expected

    with pytest.raises(ValueError):
        ccd.single_ccd("No existing name")


def test_ccd_str(ccd):
    actual = ccd.__str__()
    print(actual)
    expected = """name: Va_O1
--------------------------------------------------
name: from_0_to_1
omega: 0.10
carriers: h e
  dQ    disp ratio    corr. energy    relative energy  used for fitting?    is shallow?
0.00          0.00            3.40               3.40  True                 False
1.00          0.50            2.20               2.20  True                 False
2.00          1.00            3.30               3.30  True                 False
3.00          1.50            3.50               3.50  False                False
--------------------------------------------------
name: from_1_to_0
omega: N.A.
  dQ    disp ratio    corr. energy    relative energy  used for fitting?    is shallow?
1.00          0.80           10.30              10.30  -                    False
2.00          0.90           10.20              10.20  -                    False
3.00          1.00           10.10              10.10  -                    False"""
    assert actual == expected


def test_spline3():
    x = [0.0, 1.0, 2.0, 3.0]
    y = [0.0, 1.0, 8.0, 28.0]
    assert (spline3(x, y, num_points=11))[1][1] == pytest.approx(2.17924820)
    # actual = spline3(x, y, num_points=11, x_min=-1.0, x_max=4.0)
    # print(actual)
    # assert actual[0][0] == pytest.approx(-1.0)
    # assert actual[0][-1] == pytest.approx(4.0)


def test_plot_ccd(ccd):
    plotter = CcdPlotter(ccd)
    plotter.construct_plot()
    plotter.plt.show()
