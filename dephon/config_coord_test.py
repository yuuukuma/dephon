# -*- coding: utf-8 -*-
#  Copyright (c) 2022 Kumagai group.
import numpy as np
import pytest

from dephon.config_coord import Ccd, SinglePointInfo, CcdPlotter, \
    SingleCcd, spline3
from dephon.enum import CorrectionEnergyType, Carrier

band_edges = dict(vbm=1.0, cbm=2.0, supercell_vbm=1.1, supercell_cbm=1.9)


@pytest.fixture
def single_point_min_info():
    return SinglePointInfo(dQ=1.0, disp_ratio=0.1)


@pytest.fixture
def single_point_max_info():
    return SinglePointInfo(dQ=1.0,
                           disp_ratio=0.1,
                           energy=2.0,
                           is_shallow=True,
                           correction_energy=3.0,
                           correction_method=CorrectionEnergyType.extended_FNV,
                           used_for_fitting=True,
                           base_energy=1.0)


def test_single_point_info_corrected_energy(single_point_min_info,
                                            single_point_max_info):
    assert single_point_min_info.relative_corrected_energy is None
    assert single_point_max_info.relative_corrected_energy == 2.0 + 3.0 - 1.0


def test_image_structure_str(single_point_min_info,
                             single_point_max_info):
    actual = single_point_min_info.__str__()
    expected = """  dQ    disp ratio  energy    corr. energy    used for fitting?    is shallow?
1.00          0.10  -         -               -                    -"""
    assert actual == expected

    actual = single_point_max_info.__str__()
    expected = """  dQ    disp ratio    energy    corr. energy  used for fitting?    is shallow?
1.00          0.10      1.00            3.00  True                 True"""
    assert actual == expected


@pytest.fixture
def single_ccd():
    return SingleCcd(
        name="q=0",
        carriers=[Carrier.hole, Carrier.electron],
        point_infos=[
            SinglePointInfo(3., 1.0, 2.3, False, 1., used_for_fitting=True),
            SinglePointInfo(2., 0.9, 0.2, False, 2., used_for_fitting=True),
            SinglePointInfo(1., 0.8, 0.3, False, 3., used_for_fitting=True),
            SinglePointInfo(4., 0.8, 0.3, False, 3., used_for_fitting=False)])


def test_single_ccd_sort_single_point_infos(single_ccd):
    actual = single_ccd.point_infos[0]
    expected = SinglePointInfo(1., 0.8, 0.3, False, 3.,
                               used_for_fitting=True, base_energy=0.0)
    assert actual == expected


def test_single_ccd_dQs_and_energies(single_ccd):
    single_ccd.set_base_energy(0.1)

    actual = single_ccd.dQs_and_energies(False)
    dQs = [1.0, 2.0, 3.0, 4.0]
    energies = [3.2, 2.1, 3.2, 3.2]
    expected = (dQs, energies)
    assert np.array(actual) == pytest.approx(np.array(expected))

    single_ccd.set_base_energy()

    actual = single_ccd.dQs_and_energies(only_used_for_fitting=True)
    dQs = [1.0, 2.0, 3.0]
    energies = [1.1, 0.0, 1.1]
    expected = (dQs, energies)
    assert np.array(actual) == pytest.approx(np.array(expected))


def test_image_structure_infos_omega(single_ccd):
    assert single_ccd.omega() == pytest.approx(0.09589760387182514)


def test_image_structure_infos_str(single_ccd):
    single_ccd.set_base_energy()
    actual = single_ccd.__str__()
    print(actual)
    expected = """name: q=0
omega: 0.10
carriers: h e
  dQ    disp ratio    energy    corr. energy  used for fitting?    is shallow?
1.00          0.80     -1.90            3.00  True                 False
2.00          0.90     -2.00            2.00  True                 False
3.00          1.00      0.10            1.00  True                 False
4.00          0.80     -1.90            3.00  False                False"""
    assert actual == expected


@pytest.fixture
def ccd(single_ccd):
    excited_state = SingleCcd(
        name="q=1",
        point_infos=[SinglePointInfo(3., 1.0, 10.1, False, 1.),
                     SinglePointInfo(2., 0.9, 10.2, False, 2.),
                     SinglePointInfo(1., 0.8, 10.3, False, 3.)])
    return Ccd(defect_name="Va_O1", ccds=[single_ccd, excited_state])


def test_ccd_single_ccd(ccd, single_ccd):
    actual = ccd.single_ccd("q=0")
    expected = single_ccd
    assert actual == expected

    with pytest.raises(ValueError):
        ccd.single_ccd("No existing name")


def test_ccd_lowest_energy(ccd):
    assert ccd.lowest_energy == 2.2


def test_ccd_str(ccd):
    actual = ccd.__str__()
    print(actual)
    expected = """name: Va_O1
--------------------------------------------------
name: q=0
omega: 0.10
carriers: h e
  dQ    disp ratio    energy    corr. energy  used for fitting?    is shallow?
1.00          0.80      0.30            3.00  True                 False
2.00          0.90      0.20            2.00  True                 False
3.00          1.00      2.30            1.00  True                 False
4.00          0.80      0.30            3.00  False                False
--------------------------------------------------
name: q=1
omega: N.A.
  dQ    disp ratio    energy    corr. energy  used for fitting?    is shallow?
1.00          0.80     10.30            3.00  -                    False
2.00          0.90     10.20            2.00  -                    False
3.00          1.00     10.10            1.00  -                    False"""
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
