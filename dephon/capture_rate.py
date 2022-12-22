# -*- coding: utf-8 -*-
#  Copyright (c) 2022 Kumagai group.
from dataclasses import dataclass
from typing import Union, List

import numpy as np
from nonrad import get_C

from dephon.config_coord import Ccd
from dephon.dephon_init import DephonInit
from dephon.ele_phon_coupling import EPCoupling


@dataclass
class CaptureRate:
    Wif: float
    phonon_overlap: float
    degeneracy: int
    volume: float

    @property
    def capture_rate(self):
        return (2 * np.pi * self.degeneracy * self.Wif**2 * self.volume
                * self.phonon_overlap)


def make_capture_rate(ccd_init: DephonInit,
                      ccd: Ccd,
                      ground_name: str,
                      excited_name: str,
                      ep_coupling: EPCoupling,
                      T: List[float]):
    return get_C(dQ=ccd_init.dQ,
                 dE=ccd_init.dE,
                 wi=ccd.omega(ground_name),
                 wf=ccd.omega(excited_name),
                 Wif=ep_coupling.averaged_wif,
                 volume=ep_coupling.volume,
                 g=ccd_init.multiplicity_by_symmetry_reduction,
                 T=np.array(T))


def get_R(dQ: float,
          dE: float,
          wi: float,
          wf: float,
          T: Union[float, np.ndarray],
          occ_tol=1e-4):
    return get_C(dQ, dE, wi, wf, Wif=1.0, volume=1.0, g=1, T=T, occ_tol=occ_tol)