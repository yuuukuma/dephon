# -*- coding: utf-8 -*-
#  Copyright (c) 2022 Kumagai group.
from dataclasses import dataclass, field
from typing import List, Dict

import matplotlib
import numpy as np
from matplotlib import pyplot as plt
from matplotlib.axes import Axes
from monty.json import MSONable
from nonrad.ccd import get_omega_from_PES
from pymatgen.analysis.defects.ccd import get_dQ
from pymatgen.core import Structure
from scipy import interpolate
from tabulate import tabulate
from vise.util.logger import get_logger
from vise.util.matplotlib import float_to_int_formatter
from vise.util.mix_in import ToJsonFileMixIn

logger = get_logger(__name__)


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
class MinimumPointInfo(MSONable):
    charge: int
    structure: Structure
    energy: float
    energy_correction: float
    initial_site_symm: str
    final_site_symm: str
    site_symmetry_opt_num: int

    @property
    def corrected_energy(self):
        return self.energy + self.energy_correction


@dataclass
class CcdInit(MSONable, ToJsonFileMixIn):
    name: str
    ground_state: MinimumPointInfo
    excited_state: MinimumPointInfo
    vbm: float
    cbm: float
    supercell_vbm: float
    supercell_cbm: float
    multiplicity_by_symmetry_reduction: int = None

    def __post_init__(self):
        if abs(self.ground_state.charge - self.excited_state.charge) != 1:
            logger.warning("The charge difference is not 1. You must know what"
                           "you are doing.")

    @property
    def ground_structure(self):
        return self.ground_state.structure

    @property
    def excited_structure(self):
        return self.excited_state.structure

    @property
    def ground_charge(self):
        return self.ground_state.charge

    @property
    def excited_charge(self):
        return self.excited_state.charge

    @property
    def ground_corrected_energy(self):
        return self.ground_state.corrected_energy

    @property
    def excited_corrected_energy(self):
        return self.excited_state.corrected_energy

    @property
    def dQ(self):
        return get_dQ(self.excited_structure, self.ground_structure)

    @property
    def zero_phonon_line(self):
        return self.excited_corrected_energy - self.ground_corrected_energy

    @property
    def dR(self):
        return get_dR(self.excited_structure, self.ground_structure)

    @property
    def modal_mass(self):
        return (self.dQ / self.dR) ** 2

    @property
    def charge_change(self):
        return self.excited_charge - self.ground_charge

    def __str__(self):
        result = [f"name: {self.name}"]
        table = [["vbm", self.vbm, "supercell vbm", self.supercell_vbm],
                 ["cbm", self.cbm, "supercell cbm", self.supercell_cbm],
                 ["dQ (amu^0.5 Å)", self.dQ],
                 ["dR (Å)", self.dR],
                 ["M (amu)", self.modal_mass]]
        result.append(tabulate(table, tablefmt="plain", floatfmt=".3f"))

        result.append("-" * 60)

        headers = ["state", "initial symm", "final symm", "energy",
                   "correction", "corrected energy"]
        table = []
        for s in [self.excited_state, self.ground_state]:
            state_str = f"{self.name}_{s.charge}"
            table.append([state_str, s.initial_site_symm, s.final_site_symm,
                          s.energy, s.energy_correction, s.corrected_energy])

        result.append(
            tabulate(table, tablefmt="plain", headers=headers, floatfmt=".3f",
                     stralign="center"))
        result.append(f"ZPL: {self.zero_phonon_line:.3f}")
        return "\n".join(result)


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
    fitting_q_ranges: Dict[str, List[float]] = field(default_factory=dict)

    def __post_init__(self):
        for v in self.image_infos.values():
            v.sort(key=lambda x: x.dQ)

    @property
    def lowest_energy(self):
        energies = []
        for imag_infos in self.image_infos.values():
            energies += [imag_info.corrected_energy for imag_info in imag_infos]
        return min(energies)

    def set_q_range(self, image_name, min_q, max_q):
        if min_q is None and max_q is None:
            logger.info(f"Quadratic fitting range for {image_name} is unset.")
            self.fitting_q_ranges.pop(image_name)
            return
        if min_q is None or max_q is None:
            logger.warning(f"To set quadratic fitting range, both min and max"
                           f"values need to be set.")
            return

        logger.info(f"Quadratic fitting for {image_name} is set to "
                    f"[{min_q}, {max_q}].")
        self.fitting_q_ranges[image_name] = [min_q, max_q]

    def q_range(self, image_name):
        return self.fitting_q_ranges.get(image_name, None)

    def omega(self,
              image_name: str,
              ax: Axes = None,
              energy_std: float = 0.0):
        dQs, energies = [], []
        q_range = self.q_range(image_name)
        if q_range is None:
            raise ValueError("To calculate omega, set Q range.")

        for v in self.image_infos[image_name]:
            if q_range[0] < v.dQ < q_range[1]:
                dQs.append(v.dQ)
                energies.append(v.corrected_energy - energy_std)

        if len(dQs) < 2:
            raise ValueError("The number of Q points is not sufficient.")

        return get_omega_from_PES(np.array(dQs), np.array(energies), ax=ax)

    def __str__(self):
        result = [f"name: {self.name}"]
        headers = ["dQ", "energy", "corr", "corr energy (rel)", "corr method"]
        for image_name, imag_structures in self.image_infos.items():
            result.append(image_name)
            table = []
            for imag_structure in imag_structures:
                table.append([imag_structure.dQ,
                              imag_structure.energy,
                              imag_structure.correction,
                              imag_structure.corrected_energy - self.lowest_energy,
                              imag_structure.correction_type])
            tabu = tabulate(table, tablefmt="plain", floatfmt=".2f",
                            headers=headers)
            result.append(tabu)

            q_range = self.q_range(image_name)
            if q_range:
                q_range = f"[{q_range[0]}, {q_range[1]}]"
                result.append(f"Fitting Q range {q_range}")
                result.append(f"omega: {self.omega(image_name):.3f}")
            result.append("")

        return "\n".join(result)


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
        min_e = self._ccd.lowest_energy
        for more_name, v in self._ccd.image_infos.items():
            dQs = [vv.dQ for vv in v]
            energies = [vv.corrected_energy - min_e for vv in v]
            x, y = spline3(dQs, energies, 100, self._spline_deg)
            self.plt.plot(x, y, label=more_name)
            self.plt.scatter(dQs, energies, marker='o')

            ax = self.plt.gca()
            if self._ccd.fitting_q_ranges.get(more_name, None):
                self._ccd.omega(more_name, ax=ax, energy_std=min_e)

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

