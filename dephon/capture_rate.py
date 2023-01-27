# -*- coding: utf-8 -*-
#  Copyright (c) 2022 Kumagai group.
from dataclasses import dataclass
from typing import List

import numpy as np
from monty.json import MSONable
from nonrad import get_C
from tabulate import tabulate
from vise.util.mix_in import ToJsonFileMixIn

from dephon.config_coord import SingleCcd


@dataclass
class CaptureRate(MSONable, ToJsonFileMixIn):
    Wif: float
    phonon_overlaps: List[float]  # as a function of temperature
    temperatures: List[float]
    degeneracy: float
    volume: float

    @property
    def capture_rate(self) -> np.array:
        return (2 * np.pi * self.degeneracy * self.Wif ** 2 * self.volume
                * np.array(self.phonon_overlaps))

    def __str__(self):
        result = []
        aaa = [["Wif:", self.Wif],
               ["degeneracy:", self.degeneracy],
               ["volume (Å):", self.volume]]

        result.append(tabulate(aaa, tablefmt="plain", floatfmt=".3f"))
        result.append("phonon overlap:")

        phonon_overlaps = []
        for omega, T in zip(self.phonon_overlaps, self.temperatures):
            phonon_overlaps.append([T, omega])
        result.append(tabulate(phonon_overlaps,
                               headers=["K", "THz"],
                               tablefmt="plain", floatfmt=".3f"))
        result.append("-" * 50)
        result.append(f"capture rate:")

        cap_rates = []
        for rate, T in zip(self.capture_rate, self.temperatures):
            cap_rates.append([T, rate])
        result.append(tabulate(phonon_overlaps,
                               headers=["K", "XXX"],
                               tablefmt="plain", floatfmt=".3f"))
        return "\n".join(result)


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

