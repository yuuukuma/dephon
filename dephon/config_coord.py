# -*- coding: utf-8 -*-
#  Copyright (c) 2022 Kumagai group.
from dataclasses import dataclass
from typing import List, Dict

import numpy as np
from matplotlib import pyplot as plt
from monty.json import MSONable
from pymatgen.analysis.defects.ccd import get_dQ
from pymatgen.core import Structure
from scipy import interpolate
from tabulate import tabulate
from vise.util.matplotlib import float_to_int_formatter
from vise.util.mix_in import ToJsonFileMixIn


def get_dR(ground: Structure, excited: Structure) -> float:
    """Summation of atomic displacement distance

    Args:
        ground (Structure): Reference structure
        excited (Structure): Target structure

    Note that constant deviation is removed by e.g., aligning the farthest atom.

    Returns:
        The Summed atomic displacement distance in float
    """
    return np.sqrt(
        np.sum([x.distance(y) ** 2 for x, y in zip(ground, excited)]))


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
    vbm: float
    cbm: float
    supercell_vbm: float
    supercell_cbm: float

    @property
    def dQ(self):
        return get_dQ(self.excited_structure, self.ground_structure)

    @property
    def dR(self):
        return get_dR(self.excited_structure, self.ground_structure)

    @property
    def modal_mass(self):
        return (self.dQ / self.dR) ** 2

    def __str__(self):
        charge_diff = self.ground_charge - self.excited_charge
        if charge_diff == -1:
            trapped_carrier = "e-"
        elif charge_diff == 1:
            trapped_carrier = "h+"
        else:
            raise ValueError("The charge difference between excited and ground "
                             "states is neither 1 nor -1")

        result = [["name", self.name],
                  ["vbm", round(self.vbm, 2),
                   "supercell vbm", round(self.supercell_vbm, 2)],
                  ["cbm", round(self.cbm, 2),
                   "supercell cbm", round(self.supercell_cbm, 2)],
                  ["dQ (amu^0.5 Å)", round(self.dQ, 2)],
                  ["dR (Å)", round(self.dR, 2)],
                  ["M (amu)", round(self.modal_mass, 2)]]
        excited_state = f"{self.name}_{self.excited_charge} + {trapped_carrier}"
        result.append(["Excited state:", excited_state,
                       "energy:", self.excited_energy,
                       "correction:", round(self.excited_energy_correction, 3)])

        ground_state = f"{self.name}_{self.ground_charge}"
        result.append(["Ground state:", ground_state,
                       "energy:", self.ground_energy,
                       "correction:", round(self.ground_energy_correction, 3)])

        return tabulate(result, tablefmt="plain")


@dataclass
class ImageStructureInfo(MSONable, ToJsonFileMixIn):
    dQ: float
    energy: float = None
    correction: float = None
    correction_type: str = None

    @property
    def corrected_energy(self):
        return self.energy + self.correction


@dataclass
class Ccd(MSONable, ToJsonFileMixIn):
    image_infos: Dict[str, List[ImageStructureInfo]]
    name: str = None

    def __post_init__(self):
        for v in self.image_infos.values():
            v.sort(key=lambda x: x.dQ)

    @property
    def min_energy(self):
        energies = []
        for imag_infos in self.image_infos.values():
            energies += [imag_info.corrected_energy for imag_info in imag_infos]
        return min(energies)


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
        min_e = self._ccd.min_energy
        for more_name, v in self._ccd.image_infos.items():
            dQs = [vv.dQ for vv in v]
            energies = [vv.corrected_energy - min_e for vv in v]
            x, y = spline3(dQs, energies, 100, self._spline_deg)
            self.plt.plot(x, y, label=more_name)
            self.plt.scatter(dQs, energies, marker='o')

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

