# -*- coding: utf-8 -*-
#  Copyright (c) 2022 Kumagai group.
from copy import deepcopy
from dataclasses import dataclass, field
from math import isclose
from typing import List, Optional

import numpy as np
from matplotlib import pyplot as plt
from matplotlib.axes import Axes
from monty.json import MSONable
from nonrad.ccd import get_omega_from_PES
from nonrad.nonrad import AMU2KG, ANGS2M, EV2J, HBAR
from scipy import interpolate
from scipy.optimize import curve_fit
from tabulate import tabulate
from vise.util.logger import get_logger
from vise.util.matplotlib import float_to_int_formatter
from vise.util.mix_in import ToJsonFileMixIn

from dephon.enum import CorrectionType, Carrier

logger = get_logger(__name__)


_imag_headers = ["dQ", "disp ratio", "corr. energy", "relative energy",
                 "used for fitting?", "is shallow?"]


def get_omega_from_PES(
        Q: np.ndarray,
        energy: np.ndarray,
        Q0: Optional[float] = None,
        ax=None,
        q: Optional[np.ndarray] = None
) -> float:
    """Calculate the harmonic phonon frequency for the given PES.

    Parameters
    ----------
    Q : np.array(float)
        array of Q values (amu^{1/2} Angstrom) corresponding to each vasprun
    energy : np.array(float)
        array of energies (eV) corresponding to each vasprun
    Q0 : float
        fix the minimum of the parabola (default is None)
    ax : matplotlib.axes.Axes
        optional axis object to plot the resulting fit (default is None)
    q : np.array(float)
        array of Q values to evaluate the fitting function at

    Returns
    -------
    float
        harmonic phonon frequency from the PES in eV
    """
    def f(Q, omega, Q0, dE):
        return 0.5 * omega**2 * (Q - Q0)**2 + dE

    # set bounds to restrict Q0 to the given Q0 value
    bounds = (-np.inf, np.inf) if Q0 is None else \
        ([-np.inf, Q0 - 1e-10, -np.inf], [np.inf, Q0, np.inf])
    popt, _ = curve_fit(f, Q, energy, bounds=bounds)    # pylint: disable=W0632

    # optional plotting to check fit
    if ax is not None:
        q_L = np.max(Q) - np.min(Q)
        if q is None:
            q = np.linspace(np.min(Q) - 0.1 * q_L, np.max(Q) + 0.1 * q_L, 1000)
        ax.plot(q, f(q, *popt))

    return HBAR * popt[0] * np.sqrt(EV2J / (ANGS2M**2 * AMU2KG))


@dataclass
class SinglePointInfo(MSONable, ToJsonFileMixIn):
    dQ: float
    disp_ratio: float
    # energy obtained directly from DFT calculations
    corrected_energy: float = None
    is_shallow: bool = None
    correction_method: CorrectionType = None
    # whether this point is used for the quadratic fitting.
    used_for_fitting: bool = None
    # This needs to be here to make table_for_plot
    base_energy: float = 0.0  # This must be same in the same (Single)Ccd class

    @property
    def relative_energy(self):
        if self.corrected_energy is None:
            return None
        return self.corrected_energy - self.base_energy

    @property
    def table_for_plot(self):
        result = [self.dQ, self.disp_ratio, self.corrected_energy,
                  self.relative_energy, self.used_for_fitting, self.is_shallow]
        return result
#        return ["-" if x is None else x for x in result]

    def __str__(self):
        return tabulate([self.table_for_plot], tablefmt="plain", floatfmt=".3f",
                        headers=_imag_headers)


@dataclass
class SingleCcd(MSONable, ToJsonFileMixIn):
    name: str
    charge: int
    point_infos: List[SinglePointInfo] = field(default_factory=list)
    carriers: List[Carrier] = field(default_factory=list)

    def __post_init__(self):
        self.point_infos.sort(key=lambda x: x.dQ)

    def set_base_energy(self, energy=None):
        if energy is None:
            energy = self.disp_point_info(0.0).corrected_energy

        for i in self.point_infos:
            i.base_energy = energy

    def set_quadratic_fitting_range(self, q_range: list = None):
        q_min, q_max = q_range if q_range else [float("-inf"), float("inf")]
        for p_info in self.point_infos:
            p_info.used_for_fitting = q_min < p_info.dQ < q_max

    def shift_energy(self, energy):
        for i in self.point_infos:
            if i.corrected_energy:
                i.corrected_energy += energy

    def disp_point_info(self, disp) -> SinglePointInfo:
        for point_info in self.point_infos:
            if isclose(point_info.disp_ratio, disp):
                return point_info
        raise ValueError()

    def dQs_and_energies(self, only_used_for_fitting=False):
        dQs, energies = [], []
        for imag in self.point_infos:
            if only_used_for_fitting and imag.used_for_fitting is not True:
                continue
            if imag.relative_energy is None:
                continue
            dQs.append(imag.dQ)
            energies.append(imag.relative_energy)
        return dQs, energies

    def omega(self, ax: Axes = None):
        dQs, energies = self.dQs_and_energies(True)
        if len(dQs) < 3:
            raise ValueError("The number of Q points is not sufficient for "
                             "calculating omega.")

        return get_omega_from_PES(np.array(dQs), np.array(energies), ax=ax)

    def dQ_reverted_single_ccd(self) -> "SingleCcd":
        result = deepcopy(self)
        disp1_dQ = self.disp_point_info(1.0).dQ
        for i in result.point_infos:
            i.dQ = (1.0 - i.disp_ratio) * disp1_dQ
        result.point_infos.sort(key=lambda x: x.dQ)
        return result

    def add_plot(self, ax, color, q_range):
        dQs, energies = self.dQs_and_energies()
        ax.scatter(dQs, energies, marker='o', color=color)
        try:
            print(dQs, energies)
            x, y = spline3(dQs, energies, 100, q_range)
            ax.plot(x, y, label=self.name, color=color)
        except TypeError as e:
            print(f"{self.name}: {e}")
            pass

        try:
            self.omega(ax)
        except (ValueError, TypeError) as e:
            print(f"{self.name}: {e}")
            pass

    def __str__(self):
        try:
            omega = f"{self.omega():.3f}"
        except ValueError:
            omega = "N.A."
        result = [f"name: {self.name}", f"charge: {self.charge}",
                  f"omega: {omega}"]

        if self.carriers:
            carriers = " ".join([str(carrier) for carrier in self.carriers])
            result.append(f"carriers: {carriers}")

        table_data = [x.table_for_plot for x in self.point_infos]
        result.append(tabulate(table_data, tablefmt="plain", floatfmt=".3f",
                               headers=_imag_headers + ["omega"]))
        return "\n".join(result)


@dataclass
class Ccd(MSONable, ToJsonFileMixIn):
    defect_name: str
    ccds: List[SingleCcd]

    def single_ccd(self, ccd_name) -> SingleCcd:
        names = []
        for i in self.ccds:
            if i.name == ccd_name:
                return i
            names.append(i.name)
        raise ValueError(f"Choose state name from {' '.join(names)}")

    def __str__(self):
        result = [f"name: {self.defect_name}"]

        for imag in self.ccds:
            result.append("-"*50)
            result.append(imag.__str__())
        return "\n".join(result)


def spline3(xs, ys, num_points, xrange=None):
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

    if xrange:
        x_dist = max(xs) - min(xs)
        _min = (xrange[0] - min(xs)) / x_dist
        _max = (xrange[1] - min(xs)) / x_dist
    else:
        _min, _max = 0.0, 1.0

    print(xs, xrange, _min, _max)
    u = np.linspace(_min, _max, num=num_points, endpoint=True)
    spline = interpolate.splev(u, tck)
    return spline[0], spline[1]


class CcdPlotter:
    def __init__(self, ccd: Ccd,
                 title: str = None,
                 set_energy_zero: bool = True,
                 q_range: list = None):
        self._title = title or ""
        self._ccd = ccd
        self._set_energy_zero = set_energy_zero
        self.plt = plt
        self._q_range = q_range

    def construct_plot(self):
        self._add_ccd()
        self._set_title()
        self._set_formatter()
        self._set_labels()
        self.plt.tight_layout()

    def _add_ccd(self):
        ax = self.plt.gca()
        if self._q_range:
            ax.set_xlim(self._q_range[0], self._q_range[1])
        for imag_infos in self._ccd.ccds:
            imag_infos.add_plot(ax, "black", self._q_range)

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

