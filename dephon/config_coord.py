# -*- coding: utf-8 -*-
#  Copyright (c) 2022 Kumagai group.
from dataclasses import dataclass
from typing import List

import numpy as np
from matplotlib import pyplot as plt
from monty.json import MSONable
from pydefect.analyzer.defect_energy import DefectEnergy
from pymatgen.analysis.defects.ccd import get_dQ
from pymatgen.core import Structure
from scipy import interpolate
from tabulate import tabulate
from vise.util.matplotlib import float_to_int_formatter
from vise.util.mix_in import ToJsonFileMixIn


@dataclass
class CcdInit(MSONable, ToJsonFileMixIn):
    name: str
    excited_structure: Structure  # excited excited structure
    ground_structure: Structure  # final ground structure
    excited_charge: int
    ground_charge: int
    excited_energy: float
    excited_energy_correction: float
    ground_energy: float
    ground_energy_correction: float

    @property
    def dQ(self):
        return get_dQ(self.excited_structure, self.ground_structure)

    def __str__(self):
        charge_diff = self.ground_charge - self.excited_charge
        if charge_diff == -1:
            trapped_carrier = "e-"
        elif charge_diff == 1:
            trapped_carrier = "h+"
        else:
            raise ValueError("The charge difference between excited and ground "
                             "states is neither 1 nor -1")

        result = []
        excited_state = f"{self.name}_{self.excited_charge} + {trapped_carrier}"
        result.append(["Excited state:", excited_state,
                       "energy:", self.excited_energy,
                       "correction:", self.excited_energy_correction])

        ground_state = f"{self.name}_{self.ground_charge}"
        result.append(["Ground state:", ground_state,
                       "energy:", self.ground_energy,
                       "correction:", self.ground_energy_correction])

        return tabulate(result, tablefmt="plain")


@dataclass
class ImageStructureInfo(MSONable):
    displace_ratio: float  # 0.0: original, 1.0: counter structure
    energy: float  # include the corrections
    dQ: float


@dataclass
class Ccd(MSONable, ToJsonFileMixIn):
    dQ: float
    excited_image_infos: List[ImageStructureInfo]
    ground_image_infos: List[ImageStructureInfo]
    correction_type: str = None

    @property
    def ground_dQs(self):
        return [i.dQ for i in self.ground_image_infos]

    @property
    def ground_energies(self):
        return [i.energy for i in self.ground_image_infos]

    @property
    def excited_dQs(self):
        return [i.dQ for i in self.excited_image_infos]

    @property
    def excited_energies(self):
        return [i.energy for i in self.excited_image_infos]

    # @property
    # def excited_dQs(self):
    #     return [self.dQ * (1 - i.displace_ratio)
    #             for i in self.excited_image_infos]

def ccd_plt(ccd: Ccd):
    plt.plot(ccd.ground_dQs, ccd.ground_energies)
    plt.plot(ccd.excited_dQs, ccd.excited_energies)
    return plt


def spline3(x, y, point, deg):
    tck, u = interpolate.splprep([x, y], k=deg, s=0)
    u = np.linspace(0, 1, num=point, endpoint=True)
    spline = interpolate.splev(u, tck)
    return spline[0], spline[1]


class CcdPlotter:
    def __init__(self, ccd: Ccd,
                 title: str = None,
                 set_energy_zero: bool = True,
                 spline_deg: int = 3):
        self._title = title or ""
        self._ccd = ccd
        self._set_energy_zero = set_energy_zero
        self._spline_deg = spline_deg
        self.plt = plt

    def construct_plot(self):
        self._add_ccd()
        self._set_title()
        self._set_formatter()
        self._set_labels()
        self.plt.tight_layout()

    def _add_ccd(self):
        g_dQs, g_energies = self._ccd.ground_dQs, self._ccd.ground_energies
        e_dQs, e_energies = self._ccd.excited_dQs, self._ccd.excited_energies
        if self._set_energy_zero:
            base_energy = min(g_energies + e_energies)
            g_energies = [e - base_energy for e in g_energies]
            e_energies = [e - base_energy for e in e_energies]

        x, y = spline3(g_dQs, g_energies, 100, self._spline_deg)
        self.plt.plot(x, y, label="splprep")
        self.plt.scatter(g_dQs, g_energies, marker='o')

        x, y = spline3(e_dQs, e_energies, 100, self._spline_deg)
        self.plt.plot(x, y, label="splprep")
        self.plt.scatter(e_dQs, e_energies, marker='o')

    def _set_labels(self):
        ax = self.plt.gca()
        ax.set_xlabel("Q (amu$^{1/2}$ Å)")
        ax.set_ylabel("Energy (eV)")

    def _set_title(self):
        self.plt.gca().set_title(self._title)

    def _set_formatter(self):
        self.plt.gca().xaxis.set_major_formatter(float_to_int_formatter)
        self.plt.gca().yaxis.set_major_formatter(float_to_int_formatter)

