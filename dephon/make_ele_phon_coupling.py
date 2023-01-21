# -*- coding: utf-8 -*-
#  Copyright (c) 2022 Kumagai group.
from typing import Dict, Tuple, Optional

import numpy as np
from pymatgen.electronic_structure.core import Spin

from dephon.config_coord import Ccd, SingleCcd
from dephon.dephon_init import DephonInit
from dephon.ele_phon_coupling import EPCoupling, DefectBandId, EPMatrixElement, \
    InnerProduct
from dephon.enum import Carrier
from dephon.util import spin_to_idx


class MakeInitialEPCoupling:
    def __init__(self,
                 dephon_init: DephonInit,
                 ccd: Ccd,
                 captured_carrier: Carrier,
                 charge_for_e_p_coupling: int = None):
        self.dephon_init = dephon_init
        self.captured_carrier = captured_carrier
        self.i_single_ccd, self.f_single_ccd = \
            ccd.initial_and_final_ccd_from_captured_carrier(captured_carrier)

        ccd = self._single_ccd_for_e_p_coupling(charge_for_e_p_coupling)
        self.charge = ccd.charge
        self._ground_point = ccd.disp_point_info(0.0)
        self._min_point = self.dephon_init.min_info_from_charge(self.charge)

    def _single_ccd_for_e_p_coupling(self, charge_for_e_p_coupling
                                     ) -> SingleCcd:
        if charge_for_e_p_coupling:
            if charge_for_e_p_coupling == self.i_single_ccd.charge:
                return self.i_single_ccd
            elif charge_for_e_p_coupling == self.f_single_ccd.charge:
                return self.f_single_ccd
            else:
                raise ValueError

        return self._default_single_ccd_for_e_p_coupling()

    def _default_single_ccd_for_e_p_coupling(self):
        if abs(self.i_single_ccd.charge) < abs(self.f_single_ccd.charge):
            return self.i_single_ccd
        else:
            return self.f_single_ccd

    def make(self):
        return EPCoupling(
            charge=self.charge,
            captured_carrier=self.captured_carrier,
            volume=self.dephon_init.volume,
            ave_captured_carrier_mass=self.dephon_init.ave_hole_mass,
            ave_static_diele_const=self.dephon_init.ave_static_diele_const,
            e_p_matrix_elements=self._e_p_matrix_elements)

    @property
    def _e_p_matrix_elements(self):
        result = {}
        if self._ground_point.localized_orbitals:
            for los_by_spin, spin in zip(self._ground_point.localized_orbitals,
                                         [Spin.up, Spin.down]):
                for loc_orb in los_by_spin:
                    if self.captured_carrier.is_occupied(loc_orb.occupation):
                        continue

                    defect_band_id = DefectBandId(loc_orb.band_idx, spin)
                    result[defect_band_id] = self._matrix_elements(loc_orb, spin)

        return result

    def _matrix_elements(self, lo, spin):
        result = {}
        near_edge_states = \
            self._min_point.near_edge_states(self.captured_carrier, spin)

        for state in near_edge_states:
            eigenvalue_diff = abs(state.eigenvalue - lo.ave_energy)
            result[state.band_index] = EPMatrixElement(eigenvalue_diff,
                                                       state.kpt_index,
                                                       state.kpt_coord)

        return result


wswq_type = Dict[Optional[Tuple[int, int]], Dict[Tuple[int, int], complex]]


def add_inner_products(e_p_coupling: EPCoupling,
                       wswqs: wswq_type,
                       dQ: float,
                       used_for_fitting: bool = True):
    """    Returns
    -------
    dict(dict)
        a dict of dicts that takes keys (spin, kpoint) and (initial, final) as
        indices and maps it to a complex number"""
    for id_, val in e_p_coupling.e_p_matrix_elements.items():
        for edge_idx, ep_matrix_elem in val.items():
            spin_idx = spin_to_idx(id_.spin_channel)
            spin_kpt_pair = (spin_idx, ep_matrix_elem.kpt_idx)
            band_indices = [id_.defect_band_index, edge_idx]
            band_indices.sort()
            braket = np.abs(wswqs[spin_kpt_pair][tuple(band_indices)])
            inner_prod = InnerProduct(inner_product=braket,
                                      dQ=dQ,
                                      used_for_fitting=used_for_fitting)
            ep_matrix_elem.inner_products.append(inner_prod)
