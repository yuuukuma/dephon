# -*- coding: utf-8 -*-
#  Copyright (c) 2022 Kumagai group.

import pytest

from dephon.config_coord import Ccd, SinglePointInfo, CcdPlotter, \
    SingleCcd
from dephon.enum import Carrier


@pytest.fixture
def single_ccd():
    return SingleCcd(
        name="from_0_to_1",
        charge=0,
        carriers=[Carrier.h, Carrier.e],
        point_infos=[
            SinglePointInfo(2., 1.0, 3.3, False, used_for_fitting=True),
            SinglePointInfo(1., 0.5, 2.2, False, used_for_fitting=True),
            SinglePointInfo(0., 0.0, 3.4, False, used_for_fitting=True),
            SinglePointInfo(3., 1.5, 3.5, False, used_for_fitting=False)])


@pytest.fixture
def ccd(single_ccd):
    excited_state = SingleCcd(
        name="from_1_to_0",
        charge=1,
        point_infos=[SinglePointInfo(3., 1.0, 10.1, False),
                     SinglePointInfo(2., 0.9, 10.2, False),
                     SinglePointInfo(1., 0.8, 10.3, False)])
    return Ccd(defect_name="Va_O1", ccds=[single_ccd, excited_state])


def test_plot_ccd(ccd):
    plotter = CcdPlotter(ccd)
    plotter.construct_plot()
    plotter.plt.show()
