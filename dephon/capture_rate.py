# -*- coding: utf-8 -*-
#  Copyright (c) 2022 Kumagai group.
from dataclasses import dataclass
from typing import List

import numpy as np
from nonrad import get_C

from dephon.config_coord import Ccd, SingleCcd
from dephon.dephon_init import DephonInit
from dephon.ele_phon_coupling import EPCoupling


@dataclass
class CaptureRate:
    Wif: float
    phonon_overlap: List[float]
    degeneracy: float
    volume: float
    temperatures: List[float]

    @property
    def capture_rate(self):
        return (2 * np.pi * self.degeneracy * self.Wif**2 * self.volume
                * self.phonon_overlap)


def make_capture_rate(dephon_init: DephonInit,
                      ccd: Ccd,
                      e_p_coupling: EPCoupling,
                      temperatures: List[float]):
    carrier = e_p_coupling.captured_carrier
    i_ccd, f_ccd = ccd.initial_and_final_ccd_from_captured_carrier(carrier)
    i_min_info = dephon_init.min_info_from_charge(i_ccd.charge)
    f_min_info = dephon_init.min_info_from_charge(f_ccd.charge)
    i_deg = i_min_info.degeneracy_by_symmetry_reduction
    f_deg = f_min_info.degeneracy_by_symmetry_reduction

    return CaptureRate(e_p_coupling.wif,
                       calc_phonon_overlaps(i_ccd, f_ccd, temperatures),
                       f_deg / i_deg,
                       e_p_coupling.volume,
                       temperatures)


def calc_phonon_overlaps(ground_ccd: SingleCcd,
                         excited_ccd: SingleCcd,
                         T: List[float]):
    result = get_C(dQ=abs(excited_ccd.ground_point_info.dQ
                          - ground_ccd.ground_point_info.dQ),
                   dE=(excited_ccd.ground_point_info.corrected_energy
                       - ground_ccd.ground_point_info.corrected_energy),
                   wi=ground_ccd.omega(),
                   wf=excited_ccd.omega(),
                   Wif=1,
                   volume=1,
                   g=1,
                   T=np.array(T))
    return list(result)

