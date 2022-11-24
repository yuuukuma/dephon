# -*- coding: utf-8 -*-
#  Copyright (c) 2022 Kumagai group.
from copy import deepcopy

from pydefect.analyzer.band_edge_states import BandEdgeOrbitalInfos, OrbitalInfo

from dephon.plot_eigenvalues import EigenvaluePlotter


def test_eigenvalue_plotter_no_mag(mocker):
    i1 = mocker.Mock(spec=BandEdgeOrbitalInfos, autospec=True)
    i2 = mocker.Mock(spec=BandEdgeOrbitalInfos, autospec=True)
    i3 = mocker.Mock(spec=BandEdgeOrbitalInfos, autospec=True)
    i4 = mocker.Mock(spec=BandEdgeOrbitalInfos, autospec=True)

    i1.orbital_infos = [[[OrbitalInfo(-1.1, {}, 1.0),
                          OrbitalInfo(-1.0, {}, 1.0),
                          OrbitalInfo(0.0, {}, 1.0),
                          OrbitalInfo(1.0, {}, 0.1),
                          OrbitalInfo(2.0, {}, 0.0)]]]

    i2.orbital_infos = deepcopy(i1.orbital_infos)
    i2.orbital_infos[0][0][2].energy = 0.1

    i3.orbital_infos = deepcopy(i1.orbital_infos)
    i3.orbital_infos[0][0][2].energy = 0.2

    i4.orbital_infos = deepcopy(i1.orbital_infos)
    i4.orbital_infos[0][0][2].energy = 0.3

    plotter = EigenvaluePlotter(orb_infos=[i1, i2, i3, i4],
                                qs=[-1.0, 0.0, 1.0, 2.0])
    plotter.construct_plot()
    plotter.plt.show()


"""
TODO
1. support spin
2. raise assertion error for inappropriate calc

"""