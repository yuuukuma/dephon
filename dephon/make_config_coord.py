# -*- coding: utf-8 -*-
#  Copyright (c) 2022 Kumagai group.
from copy import deepcopy

from vise.util.logger import get_logger

from dephon.config_coord import SingleCcd, Ccd
from dephon.dephon_init import DephonInit
from dephon.enum import Carrier, BandEdge

logger = get_logger(__name__)


class MakeCcd:
    """
    Define charge_diff = q_{excited} - q_{ground}.
    - charge_diff = +1
      - Fermi level locates at the VBM (p-type)
      - carriers are recombined via (ground + h + e) → (excited + h) → ground
      - minority carrier is e and majority carrier is h

    - charge_diff = -1
      - Fermi level locates at the CBM (n-type)
      - carriers are recombined via (ground + h + e) → (excited + e) → ground
      - minority carrier is h and majority carrier is e
    """
    def __init__(self,
                 ground_ccd: SingleCcd,
                 excited_ccd: SingleCcd,
                 dephon_init: DephonInit):
        if abs(ground_ccd.charge - excited_ccd.charge) != 1:
            raise AssertionError(
                f"The charge difference needs to be 1. Now, ground state "
                f"{ground_ccd.charge} and excited state {excited_ccd.charge}.")
        self.dephon_init = dephon_init
        self.orig_ground_ccd = ground_ccd
        self.orig_excited_ccd = excited_ccd

        self.orig_ground_ccd.shift_energy(
            ground_ccd.charge * self.band_edge_level)
        self.orig_excited_ccd.shift_energy(
            excited_ccd.charge * self.band_edge_level)

        self.ref_energy = ground_ccd.disp_point_info(0.0).corrected_energy

    @property
    def charge_diff(self) -> int:
        return self.orig_excited_ccd.charge - self.orig_ground_ccd.charge

    @property
    def carrier_coexisting_with_excited(self):
        """ This should be a majority carrier."""
        return Carrier.from_carrier_charge(- self.charge_diff)

    @property
    def band_edge(self) -> BandEdge:
        if self.charge_diff == 1:
            return BandEdge.cbm
        else:
            return BandEdge.vbm

    @property
    def band_edge_level(self) -> float:
        return self.dephon_init.__getattribute__(str(self.band_edge))

    def carrier_energy(self, carrier: Carrier) -> float:
        if carrier == Carrier.electron:
            band_edge = self.dephon_init.cbm
        elif carrier == Carrier.hole:
            band_edge = self.dephon_init.vbm
        else:
            raise ValueError
        return - carrier.charge * band_edge

    @property
    def _ground_ccd(self) -> SingleCcd:
        """
        Returns: Ground state ccd with the lowest energy to be zero
        """
        result = deepcopy(self.orig_ground_ccd)
        result.set_base_energy(self.ref_energy)
        result.name = "ground"
        return result

    @staticmethod
    def _add_carrier(ccd: SingleCcd, carrier: Carrier) -> None:
        ccd.carriers.append(carrier)
        ccd.name += f" + {carrier}"

    @property
    def _ground_pn_ccd(self) -> SingleCcd:
        """
        Returns: Ground state ccd with an excited electron at CBM from VBM
        """
        result = deepcopy(self.orig_ground_ccd)
        result.set_base_energy(self.ref_energy)
        result.name = "ground"

        self._add_carrier(result, Carrier.hole)
        self._add_carrier(result, Carrier.electron)
        result.shift_energy(self.dephon_init.cbm - self.dephon_init.vbm)
        return result

    @property
    def _excited_ccd(self) -> SingleCcd:
        """
        Returns: Excited state ccd capturing a minority carrier with a majority
                 carrier. The Q values are reverted and zero is set to that of
                 ground state ccd.
        """
        result = self.orig_excited_ccd.dQ_reverted_single_ccd()
        result.set_base_energy(self.ref_energy)
        result.name = "excited"
        self._add_carrier(result, self.carrier_coexisting_with_excited)
        return result

    @property
    def ccd(self) -> Ccd:
        ccds = [self._ground_ccd, self._excited_ccd, self._ground_pn_ccd]
        for ccd in ccds:
            ccd.set_quadratic_fitting_range()

        return Ccd(defect_name=self.dephon_init.defect_name, ccds=ccds)
