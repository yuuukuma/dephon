# -*- coding: utf-8 -*-
#  Copyright (c) 2022 Kumagai group.
from dataclasses import dataclass
from typing import List

import numpy as np
from matplotlib import pyplot as plt
from matplotlib.axes import Axes
from monty.json import MSONable
from nonrad.ccd import get_omega_from_PES
from scipy import interpolate
from tabulate import tabulate
from vise.util.logger import get_logger
from vise.util.matplotlib import float_to_int_formatter
from vise.util.mix_in import ToJsonFileMixIn

from dephon.enum import CorrectionEnergyType

logger = get_logger(__name__)


_imag_headers = ["dQ", "disp ratio", "energy", "corr. energy",
                 "used for fitting?", "is shallow?"]


@dataclass
class ImageStructureInfo(MSONable, ToJsonFileMixIn):
    dQ: float
    disp_ratio: float
    energy: float = None
    correction_energy: float = None
    used_for_fitting: bool = None
    is_shallow: bool = None

    @property
    def corrected_energy(self):
        try:
            return self.energy + self.correction_energy
        except TypeError:
            return None

    @property
    def list_data(self):
        result = [self.dQ, self.disp_ratio, self.energy, self.correction_energy,
                  self.used_for_fitting, self.is_shallow]
        return ["-" if x is None else x for x in result]

    def __str__(self):
        return tabulate([self.list_data], tablefmt="plain", floatfmt=".2f",
                        headers=_imag_headers)


@dataclass
class ImageStructureInfos(MSONable, ToJsonFileMixIn):
    state_name: str
    image_structure_infos: List[ImageStructureInfo]

    def __post_init__(self):
        self.image_structure_infos.sort(key=lambda x: x.dQ)

    @property
    def lowest_energy(self):
        return min([i.corrected_energy for i in self.image_structure_infos])

    def set_q_range(self, min_q=-float("inf"), max_q=float("inf")):
        for i in self.image_structure_infos:
            i.used_for_fitting = (min_q <= i.dQ <= max_q)

    def dQs_and_energies(self, check_fitting=False, base_energy=0.0):
        dQs, energies = [], []
        for imag in self.image_structure_infos:
            if check_fitting and imag.used_for_fitting is not True:
                continue
            dQs.append(imag.dQ)
            energies.append(imag.corrected_energy - base_energy)
        return dQs, energies

    def omega(self, ax: Axes = None, base_energy=0.0):
        dQs, energies = self.dQs_and_energies(True, base_energy)
        if len(dQs) < 2:
            raise ValueError("The number of Q points is not sufficient.")

        return get_omega_from_PES(np.array(dQs), np.array(energies), ax=ax)

    def add_plot(self, ax, color, base_energy=0.0):
        dQs, energies = self.dQs_and_energies(base_energy=base_energy)
        ax.scatter(dQs, energies, marker='o', color=color)
        try:
            x, y = spline3(dQs, energies, 100)
            ax.plot(x, y, label=self.state_name, color=color)
        except TypeError:
            pass

        try:
            self.omega(ax, base_energy)
        except ValueError as e:
            print(e)
            pass

    def __str__(self):
        try:
            omega = f"{self.omega():.2f}"
        except ValueError:
            omega = "N.A."
        result = [f"state: {self.state_name}", f"omega: {omega}"]

        table_data = [x.list_data for x in self.image_structure_infos]
        result.append(tabulate(table_data, tablefmt="plain", floatfmt=".2f",
                               headers=_imag_headers + ["omega"]))
        return "\n".join(result)


@dataclass
class Ccd(MSONable, ToJsonFileMixIn):
    defect_name: str
    correction_energy_type: CorrectionEnergyType
    image_infos_list: List[ImageStructureInfos]

    @property
    def image_structure_info_by_name(self, state_name):
        names = []
        for i in self.image_infos_list:
            if i.state_name == state_name:
                return i
            names.append(i.state_name)
        raise ValueError(f"Choose state name from {' '.join(names)}")

    @property
    def lowest_energy(self):
        return min([i.lowest_energy for i in self.image_infos_list])

    def __str__(self):
        result = [f"name: {self.defect_name}",
                  f"correction energy type: {self.correction_energy_type}"]

        for imag in self.image_infos_list:
            result.append("-"*50)
            result.append(imag.__str__())
        return "\n".join(result)


def spline3(xs, ys, num_points):
    """Find the B-spline representation with 3 degree of the spline.

    Args:
        xs (array_like):
        ys (array_like):
        num_points (int): Number of interpolated points including end points.

    Returns:
        Tuple of
    """
    #   tck : tuple
    #         (t,c,k) a tuple containing the vector of knots, the B-spline
    #         coefficients, and the degree of the spline.
    tck = interpolate.splprep([xs, ys])[0]
    u = np.linspace(0, 1, num=num_points, endpoint=True)
    spline = interpolate.splev(u, tck)
    return spline[0], spline[1]


class CcdPlotter:
    def __init__(self, ccd: Ccd,
                 title: str = None,
                 set_energy_zero: bool = True):
        self._title = title or ""
        self._ccd = ccd
        self._set_energy_zero = set_energy_zero
        self.plt = plt

    def construct_plot(self):
        self._add_ccd()
        self._set_title()
        self._set_formatter()
        self._set_labels()
        self.plt.tight_layout()

    def _add_ccd(self):
        ax = self.plt.gca()
        for imag_infos in self._ccd.image_infos_list:
            imag_infos.add_plot(ax, "black", self._ccd.lowest_energy)

    def _set_labels(self):
        ax = self.plt.gca()
        ax.set_xlabel("Q (amu$^{1/2}$ Å)")
        ax.set_ylabel("Energy (eV)")
        ax.legend()

    def _set_title(self):
        self.plt.gca().set_title(self._title)

    def _set_formatter(self):
        self.plt.gca().xaxis.set_major_formatter(float_to_int_formatter)
        self.plt.gca().yaxis.set_major_formatter(float_to_int_formatter)

