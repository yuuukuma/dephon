# -*- coding: utf-8 -*-
#  Copyright (c) 2022 Kumagai group.

from vise.util.logger import get_logger

from dephon.config_coord import SingleCcd, Ccd
from dephon.dephon_init import DephonInit, MinimumPointInfo
from dephon.enum import Carrier

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
        assert abs(ground_ccd.charge - excited_ccd.charge) == 1
        self.dephon_init = dephon_init
        self.orig_ground_ccd = ground_ccd
        self.orig_excited_ccd = excited_ccd

    @property
    def ground_min_point(self) -> MinimumPointInfo:
        return self.dephon_init.min_point_info_from_charge(
            self.orig_ground_ccd.charge)

    @property
    def excited_min_point(self) -> MinimumPointInfo:
        return self.dephon_init.min_point_info_from_charge(
            self.orig_excited_ccd.charge)

    @property
    def charge_diff(self) -> int:
        return self.excited_min_point.charge - self.ground_min_point.charge

    @property
    def carrier_coexisting_with_excited(self):
        """ This should be a majority carrier."""
        return Carrier.from_carrier_charge(self.charge_diff == 1)

    def carrier_energy(self, carrier: Carrier) -> float:
        return self.dephon_init.band_gap if carrier == Carrier.electron else 0.0

    def _add_carrier(self, ccd: SingleCcd, carrier: Carrier) -> None:
        ccd.carriers.append(carrier)
        ccd.shift_energy(self.carrier_energy(carrier))
        ccd.name += f" + {carrier}"

    @property
    def _ground_ccd(self) -> SingleCcd:
        """
        Returns: Ground state ccd with the lowest energy to be zero
        """
        min_energy = self.ground_min_point.corrected_energy
        result = self.orig_ground_ccd.energy_shifted_single_ccd(min_energy)
        result.name = "ground"
        return result

    @property
    def _ground_pn_ccd(self) -> SingleCcd:
        """
        Returns: Ground state ccd with an excited electron at CBM from VBM
        """
        result = self._ground_ccd
        self._add_carrier(result, Carrier.hole)
        self._add_carrier(result, Carrier.electron)
        return result

    @property
    def _excited_ccd(self) -> SingleCcd:
        """
        Returns: Excited state ccd capturing a minority carrier with a majority
                 carrier. The Q values are reverted and zero is set to that of
                 ground state ccd.
        """
        min_energy = self.excited_min_point.corrected_energy
        result = self.orig_excited_ccd.dQ_reverted_single_ccd()
        result = result.energy_shifted_single_ccd(min_energy)
        result.name = "excited"
        self._add_carrier(result, self.carrier_coexisting_with_excited)
        return result

    @property
    def ccd(self) -> Ccd:
        return Ccd(
            defect_name=self.dephon_init.defect_name,
            ccds=[self._ground_ccd, self._excited_ccd, self._ground_pn_ccd])
