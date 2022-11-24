# -*- coding: utf-8 -*-
#  Copyright (c) 2022 Kumagai group.
from typing import List, Dict

from matplotlib import pyplot as plt
from pydefect.analyzer.band_edge_states import BandEdgeOrbitalInfos
from vise.util.logger import get_logger
from vise.util.matplotlib import float_to_int_formatter

logger = get_logger(__name__)


class EigenvaluePlotter:
    def __init__(self,
                 orb_infos: List[BandEdgeOrbitalInfos],
                 qs: List[float],
                 title: str = None):
        # try:
        #     num_kpt = list(orb_infos.values())[0].kpt_weights
        #     assert num_kpt == 1
        # except AssertionError:
        #     logger.warning(f"Currently {self.__class__} supports only single "
        #                    f"k-point.")
        #     raise

        # try:
        #     first_num_kpt_coords = orb_infos[0].kpt_coords
        #     for orb_info in orb_infos.values():
        #         assert first_num_kpt_coords == orb_info.kpt_coords
        # except AssertionError:
        #
        #     kpt_coords =

        self._title = title or ""
        self._orb_infos = orb_infos
        self._qs = qs
        self.plt = plt

    def construct_plot(self):
        self._add_eigenvalues()
        self._set_title()
        self._set_formatter()
        self._set_labels()
        self.plt.tight_layout()

    def _add_eigenvalues(self):
        num_spin = len(self._orb_infos[0].orbital_infos)
        _, axs = self.plt.subplots(num_spin)
        if num_spin == 1:
            axs = [axs]

        for spin_idx in range(num_spin):
            for orb_infos, q in zip(self._orb_infos, self._qs):
                orb_info = orb_infos.orbital_infos[spin_idx][0]
                xs = [q] * len(orb_info)
                energies = [oi.energy for oi in orb_info]
                occupations = [oi.occupation for oi in orb_info]
                axs[spin_idx].scatter(xs, energies, c=occupations,
                                      cmap='viridis')
        # self.plt.colorbar(ax=axs[0])

    def _set_labels(self):
        ax = self.plt.gca()
        ax.set_xlabel("Q (amu$^{1/2}$ Å)")
        ax.set_ylabel("Energy (eV)")

    def _set_title(self):
        self.plt.gca().set_title(self._title)

    def _set_formatter(self):
        self.plt.gca().xaxis.set_major_formatter(float_to_int_formatter)
        self.plt.gca().yaxis.set_major_formatter(float_to_int_formatter)

