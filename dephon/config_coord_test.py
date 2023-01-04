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
    expected = """   dQ    disp ratio  corr. energy    relative energy    used for fitting?    is shallow?
1.000         0.100"""
    assert actual == expected

    actual = single_point_max_info.__str__()
    expected = """   dQ    disp ratio    corr. energy    relative energy  used for fitting?    is shallow?
1.000         0.100           5.000              4.000  True                 True"""
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


def test_single_ccd_set_quadratic_fitting_range(single_ccd):
    single_ccd.set_quadratic_fitting_range(q_range=[-0.1, 0.1])
    actual = [p_info.used_for_fitting for p_info in single_ccd.point_infos]
    expected = [True, False, False, False]
    assert actual == expected

    single_ccd.set_quadratic_fitting_range()
    actual = [p_info.used_for_fitting for p_info in single_ccd.point_infos]
    expected = [True, True, True, True]
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
    actual = deepcopy(single_ccd)
    actual.shift_energy(energy=1.0)
    assert actual.point_infos[0].corrected_energy == 4.4


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
    expected = """name: from_0_to_1
charge: 0
omega: 0.098
carriers: h e
   dQ    disp ratio    corr. energy    relative energy  used for fitting?    is shallow?
0.000         0.000           3.400              0.000  True                 False
1.000         0.500           2.200             -1.200  True                 False
2.000         1.000           3.300             -0.100  True                 False
3.000         1.500           3.500              0.100  False                False"""
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
charge: 0
omega: 0.098
carriers: h e
   dQ    disp ratio    corr. energy    relative energy  used for fitting?    is shallow?
0.000         0.000           3.400              3.400  True                 False
1.000         0.500           2.200              2.200  True                 False
2.000         1.000           3.300              3.300  True                 False
3.000         1.500           3.500              3.500  False                False
--------------------------------------------------
name: from_1_to_0
charge: 1
omega: N.A.
   dQ    disp ratio    corr. energy    relative energy  used for fitting?    is shallow?
1.000         0.800          10.300             10.300                       False
2.000         0.900          10.200             10.200                       False
3.000         1.000          10.100             10.100                       False"""
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


def test_plot_ccd_q_range(ccd):
    plotter = CcdPlotter(ccd, q_range=[-2, 4])
    plotter.construct_plot()
    plotter.plt.show()
