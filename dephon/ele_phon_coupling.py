# -*- coding: utf-8 -*-
#  Copyright (c) 2022 Kumagai group.
from dataclasses import dataclass, field
from typing import List, Dict

from monty.json import MSONable
from pymatgen.electronic_structure.core import Spin
from vise.util.mix_in import ToJsonFileMixIn

from dephon.enum import Carrier


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
    eigenvalue_diff: float
    kpt_idx: int
    # Currently, symmetry is assumed not to be changed depending on dQ.
    kpt_coord: List[float]
    inner_products: List[InnerProduct] = field(default_factory=list)
    #
    # @property
    # def e_p_matrix_element(self) -> float:
    #     """ Evaluated by computing the slope of inner products"""
    #     # np.polyfit(Q, matels[i, :], 1)[0])
    #     pass


@dataclass(frozen=True)
class DefectBandId(MSONable):
    """ E-P matrix elements for _default_single_ccd_for_e_p_coupling particular defect
        band_edge_index: The band edge index starting from 1.
    """
    defect_band_index: int
    spin_channel: Spin


@dataclass
class EPCoupling(MSONable, ToJsonFileMixIn):
    """ E-P coupling constants between defect band and multiple band edges

    Attributes:
        defect_band_index: The defect band index starting from 1.
    """
    charge: int
    captured_carrier: Carrier
    volume: float
    ave_captured_carrier_mass: float
    ave_static_diele_const: float
    # int is band_edge_index.
    e_p_matrix_elements: Dict[DefectBandId, Dict[int, EPMatrixElement]] = None

    def set_inner_prod(self, d: Dict[DefectBandId, Dict[int, List[InnerProduct]]]):
        for k, v in d.items():
            for kk, vv in v.items():
                self.e_p_matrix_elements[k][kk].inner_products = vv

