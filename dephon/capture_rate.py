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
    summed_phonon_overlaps: List[float]  # as a function of temperature
    temperatures: List[float]
    site_degeneracy: float
    spin_selection_factor: float
#    velocity: float  # characteristic carrier velocity in [cm / s]
    volume: float  # in [cm3]

    @property
    def capture_rate(self) -> np.array:
        return (2 * np.pi * self.site_degeneracy * self.Wif ** 2 * self.volume
                * np.array(self.summed_phonon_overlaps))

    def __str__(self):
        aaa = [["Wif:", f"{self.Wif:.1e}"],
               ["site degeneracy:", f"{self.site_degeneracy}"],
               ["spin selection factor:", f"{self.spin_selection_factor}"],
               ["volume (Å):", f"{self.volume}"]]

        result = [tabulate(aaa, tablefmt="plain")]

        table = []
        for T, omega, rate in zip(self.temperatures, self.summed_phonon_overlaps, self.capture_rate):
            table.append([T, omega, rate])
        result.append(tabulate(table,
                               headers=["K", "THz", "capture coeff [cm3/s]"],
                               tablefmt="plain", floatfmt=[".1f", ".3f", ".1e"]))

        return "\n".join(result)


def calc_phonon_overlaps(ground_ccd: SingleCcd,
                         excited_ccd: SingleCcd,
                         T: List[float]):
    dQ = excited_ccd.ground_point_info.dQ - ground_ccd.ground_point_info.dQ

    dE = (excited_ccd.ground_point_info.corrected_energy
          - ground_ccd.ground_point_info.corrected_energy)

    result = get_C(dQ=abs(dQ),
                   dE=dE,
                   wi=ground_ccd.omega(),
                   wf=excited_ccd.omega(),
                   T=np.array(T),
                   Wif=1, volume=1, g=1)
    print(result)
    return list(result)

