# -*- coding: utf-8 -*-
#  Copyright (c) 2022 Kumagai group.
from dataclasses import dataclass, field
from typing import List, Union, Dict, Optional

import numpy as np
from monty.json import MSONable
from numpy.linalg import LinAlgError
from pymatgen.electronic_structure.core import Spin
from tabulate import tabulate
from vise.util.mix_in import ToJsonFileMixIn

from dephon.enum import Carrier
from dephon.util import spin_to_idx


@dataclass
class InnerProduct(MSONable):
    """
    Attributes:
        inner_products: \bra_{psi_i(0)} | S(0) |\ket_{psi_f(Q)}
        used_for_fitting:
    """
    inner_product: float
    used_for_fitting: bool = True


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
    # key dQ:
    inner_products: Dict[float, InnerProduct] = field(default_factory=dict)

    def __post_init__(self):
        if isinstance(self.spin, str):
            self.spin = Spin[self.spin]
        self.inner_products = \
            {float(k): v for k, v in self.inner_products.items()}

    def as_dict(self) -> dict:
        result = super().as_dict()
        result["spin"] = result["spin"].name
        return result

    @property
    def spin_idx(self):
        return spin_to_idx(self.spin)

    @property
    def dQs(self):
        return [dQ for dQ in self.inner_products.keys()]

    @property
    def inner_prods(self):
        return [ip.inner_product for ip in self.inner_products.values()]

    @property
    def _inner_prod_vs_q(self):
        return self.dQs, self.inner_prods

    def e_p_matrix_element(self, ax=None) -> Optional[float]:
        """ Evaluated by computing the slope of inner products"""
        try:
            grad, const = np.polyfit(self.dQs, self.inner_prods, 1)

            if ax:
                ax.scatter(self.dQs, self.inner_prods)

                x = np.arange(min(self.dQs), max(self.dQs), 0.01)
                y = x * grad + const
                ax.plot(x, y, alpha=0.5)

            return grad
        except (TypeError, LinAlgError):
            return None

    def __str__(self):
        result = []
        table = [["band edge index", self.band_edge_index],
                 ["defect band index", self.defect_band_index],
                 ["spin", self.spin.name],
                 ["eigenvalue difference", round(self.eigenvalue_diff, 3)],
                 ["kpoint coord", self.kpt_coord],
                 ["e-p matrix element", self.e_p_matrix_element()]]
        result.append(tabulate(table, tablefmt="plain", floatfmt=".3f"))

        inner_prods = []
        for dQ, ip in self.inner_products.items():
            inner_prods.append([dQ, ip.inner_product, ip.used_for_fitting])

        result.append(tabulate(inner_prods,
                               headers=["dQ", "inner product",
                                        "used for fitting?"],
                               tablefmt="plain", floatfmt=".3f"))

        return "\n".join(result)


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
        mass = round(self.ave_captured_carrier_mass, 2)
        diele_const = round(self.ave_static_diele_const, 2)
        table = [["charge", self.charge],
                 ["base disp", self.disp],
                 ["captured carrier", self.captured_carrier],
                 ["volume", round(self.volume, 2)],
                 ["averaged carrier mass", mass],
                 ["averaged static dielectric constant", diele_const]]
        result.append(tabulate(table, tablefmt="plain", floatfmt=".3f"))

        if self.e_p_matrix_elements is None:
            return "\n".join(result)

        for ep_elem in self.e_p_matrix_elements:
            result.append("-" * 50)
            result.append(ep_elem.__str__())

        return "\n".join(result)

    @property
    def wif(self):
        return sum([ep_elem.e_p_matrix_element()
                    for ep_elem in self.e_p_matrix_elements])
