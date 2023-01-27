# -*- coding: utf-8 -*-
#  Copyright (c) 2022 Kumagai group.
from dataclasses import dataclass, field
from typing import List, Union

import numpy as np
from monty.json import MSONable
from pymatgen.electronic_structure.core import Spin
from tabulate import tabulate
from vise.util.mix_in import ToJsonFileMixIn

from dephon.enum import Carrier


@dataclass
class Wif:
    value: float
    band_edge_index: int
    defect_band_index: int
    spin: Spin


@dataclass
class InnerProduct(MSONable):
    """
    Attributes:
        inner_products: \bra_{psi_i(0)} | S(0) |\ket_{psi_f(Q)}
        dQ:
        used_for_fitting:
    """
    inner_product: float
    dQ: float
    used_for_fitting: bool = None


@dataclass
class EPMatrixElement(MSONable):
    """ Electron-phonon (E-P) matrix element between defect band and band edge

    The phonon mode is only 1D mode.

    Attributes:
        band_edge_index: The band index for band edges starting from 1.
        eigenvalue_diff: The eigenvalue difference from final to initial states
        inner_products: List of \bra_{psi_i(0)} | S(0) |\ket_{psi_f(Q)}
            at the given Q points.
    """
    band_edge_index: int
    defect_band_index: int
    spin: Union[Spin, str]
    eigenvalue_diff: float
    kpt_idx: int
    # Currently, symmetry is assumed not to be changed depending on dQ.
    kpt_coord: List[float]
    inner_products: List[InnerProduct] = field(default_factory=list)

    def __post_init__(self):
        if isinstance(self.spin, str):
            self.spin = Spin[self.spin]

    def as_dict(self) -> dict:
        result = super().as_dict()
        result["spin"] = result["spin"].name
        return result

    @property
    def dQs(self):
        return [ip.dQ for ip in self.inner_products]

    @property
    def inner_prods(self):
        return [ip.inner_product for ip in self.inner_products]

    @property
    def _inner_prod_vs_q(self):
        return self.dQs, self.inner_prods

    def e_p_matrix_element(self, ax=None) -> float:
        """ Evaluated by computing the slope of inner products"""
        grad, const = np.polyfit(self.dQs, self.inner_prods, 1)

        if ax:
            ax.scatter(self.dQs, self.inner_prods)

            x = np.arange(min(self.dQs), max(self.dQs), 0.01)
            y = x * grad + const
            ax.plot(x, y, alpha=0.5)

        return grad


@dataclass
class EPCoupling(MSONable, ToJsonFileMixIn):
    """ E-P coupling constants between defect band and multiple band edges

    Attributes:
        defect_band_index: The defect band index starting from 1.
    """
    charge: int
    disp: float
    captured_carrier: Carrier
    volume: float
    ave_captured_carrier_mass: float
    ave_static_diele_const: float
    # int is band_edge_index.
    e_p_matrix_elements: List[EPMatrixElement] = None

    def reset_inner_prod(self):
        self.e_p_matrix_elements = []

    def __str__(self):
        result = []
        table = [["charge", self.charge],
                 ["disp", self.disp],
                 ["captured carrier", self.captured_carrier],
                 ["volume", self.volume],
                 ["averaged carrier mass", self.ave_captured_carrier_mass],
                 ["averaged static dielectric constant",
                  self.ave_static_diele_const]]
        result.append(tabulate(table, tablefmt="plain", floatfmt=".3f"))

        if self.e_p_matrix_elements is None:
            return "\n".join(result)

        result.append("-" * 60)
        header = ["band edge index",
                  "defect band index",
                  "spin",
                  "eigenvalue difference",
                  "kpoint coord",
                  "e-p matrix element"]

        e_p_table = []
        for ep_elem in self.e_p_matrix_elements:
            e_p_table.append([ep_elem.band_edge_index,
                              ep_elem.defect_band_index,
                              ep_elem.spin.name,
                              ep_elem.eigenvalue_diff,
                              ep_elem.kpt_coord,
                              ep_elem.e_p_matrix_element()])
        result.append(tabulate(e_p_table,
                               headers=header,
                               tablefmt="plain", floatfmt=".3f"))
        return "\n".join(result)

    @property
    def wif(self):
        return sum([ep_elem.e_p_matrix_element()
                    for ep_elem in self.e_p_matrix_elements])
