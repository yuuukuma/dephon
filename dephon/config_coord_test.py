# -*- coding: utf-8 -*-
#  Copyright (c) 2022 Kumagai group.

import pytest

from dephon.config_coord import Ccd, ImageStructureInfo, CcdPlotter, \
    ImageStructureInfos, spline3
from dephon.enum import CorrectionEnergyType

band_edges = dict(vbm=1.0, cbm=2.0, supercell_vbm=1.1, supercell_cbm=1.9)

@pytest.fixture
def imag_structure_info_min_info():
    return ImageStructureInfo(dQ=1.0, disp_ratio=0.1)

@pytest.fixture
def imag_structure_info_max_info():
    return ImageStructureInfo(dQ=1.0,
                              disp_ratio=0.1,
                              energy=2.0,
                              correction_energy=3.0,
                              used_for_fitting=True,
                              is_shallow=True)


def test_image_structure_info_corrected_energy(imag_structure_info_min_info,
                                               imag_structure_info_max_info):
    assert imag_structure_info_min_info.corrected_energy is None
    assert imag_structure_info_max_info.corrected_energy == 2.0 + 3.0


def test_image_structure_str(imag_structure_info_min_info,
                             imag_structure_info_max_info):
    actual = imag_structure_info_min_info.__str__()
    expected = """  dQ    disp ratio  energy    corr. energy    used for fitting?    is shallow?
1.00          0.10  -         -               -                    -"""
    assert actual == expected

    actual = imag_structure_info_max_info.__str__()
    expected = """  dQ    disp ratio    energy    corr. energy  used for fitting?    is shallow?
1.00          0.10      2.00            3.00  True                 True"""
    assert actual == expected


@pytest.fixture
def ground_state():
    return ImageStructureInfos(
        state_name="q=0",
        image_structure_infos=[ImageStructureInfo(3., 1.0, 2.3, 1.),
                               ImageStructureInfo(2., 0.9, 0.2, 2.),
                               ImageStructureInfo(1., 0.8, 0.3, 3.),
                               ImageStructureInfo(4., 0.8, 0.3, 3.),
                               ])


def test_image_structure_infos_sort(ground_state):
    actual = ground_state.image_structure_infos[0]
    expected = ImageStructureInfo(1., 0.8, 0.3, 3.)
    assert actual == expected


def test_image_structure_infos_lowest_energy(ground_state):
    assert ground_state.lowest_energy == 2.2


def test_image_structure_infos_set_q_range(ground_state):
    ground_state.set_q_range(min_q=1.1, max_q=2.1)
    actual = [i.used_for_fitting for i in ground_state.image_structure_infos]
    assert actual == [False, True, False, False]

    ground_state.set_q_range()
    actual = [i.used_for_fitting for i in ground_state.image_structure_infos]
    assert actual == [True, True, True, True]


def test_image_structure_infos_omega(ground_state):
    ground_state.set_q_range()
    assert ground_state.omega() == pytest.approx(0.04794880190203623)


def test_image_structure_infos_str(ground_state):
    ground_state.set_q_range()
    actual = ground_state.__str__()
    print(actual)
    expected = """state: q=0
omega: 0.05
  dQ    disp ratio    energy    corr. energy  used for fitting?    is shallow?
1.00          0.80      0.30            3.00  True                 -
2.00          0.90      0.20            2.00  True                 -
3.00          1.00      2.30            1.00  True                 -
4.00          0.80      0.30            3.00  True                 -"""
    assert actual == expected


@pytest.fixture
def ccd(ground_state):
    excited_state = ImageStructureInfos(
        state_name="q=1",
        image_structure_infos=[ImageStructureInfo(3., 1.0, 10.1, 1.),
                               ImageStructureInfo(2., 0.9, 10.2, 2.),
                               ImageStructureInfo(1., 0.8, 10.3, 3.)])
    return Ccd(defect_name="Va_O1",
               correction_energy_type=CorrectionEnergyType.extended_FNV,
               image_infos_list=[ground_state, excited_state])


def test_ccd_str(ccd):
    actual = ccd.__str__()
    expected = """name: Va_O1
correction energy type: extended FNV
--------------------------------------------------
state: q=0
omega: N.A.
  dQ    disp ratio    energy    corr. energy  used for fitting?    is shallow?
1.00          0.80      0.30            3.00  -                    -
2.00          0.90      0.20            2.00  -                    -
3.00          1.00      2.30            1.00  -                    -
4.00          0.80      0.30            3.00  -                    -
--------------------------------------------------
state: q=1
omega: N.A.
  dQ    disp ratio    energy    corr. energy  used for fitting?    is shallow?
1.00          0.80     10.30            3.00  -                    -
2.00          0.90     10.20            2.00  -                    -
3.00          1.00     10.10            1.00  -                    -"""
    assert actual == expected


def test_ccd_lowest_energy(ccd):
    assert ccd.lowest_energy == 2.2


def test_spline3():
    x = [0.0, 1.0, 2.0, 3.0]
    y = [0.0, 1.0, 8.0, 28.0]
    assert (spline3(x, y, num_points=11))[1][1] == pytest.approx(2.17924820)
    # actual = spline3(x, y, num_points=11, x_min=-1.0, x_max=4.0)
    # print(actual)
    # assert actual[0][0] == pytest.approx(-1.0)
    # assert actual[0][-1] == pytest.approx(4.0)


def test_plot_ccd(ccd):
    ccd.image_infos_list[0].set_q_range()
    plotter = CcdPlotter(ccd)
    plotter.construct_plot()
    plotter.plt.show()
