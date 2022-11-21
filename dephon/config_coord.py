# -*- coding: utf-8 -*-
#  Copyright (c) 2022 Kumagai group.
from dataclasses import dataclass
from typing import List

from matplotlib import pyplot as plt
from monty.json import MSONable
from pydefect.analyzer.defect_energy import DefectEnergy
from pymatgen.analysis.defects.ccd import get_dQ
from pymatgen.core import Structure
from tabulate import tabulate
from vise.util.mix_in import ToJsonFileMixIn


@dataclass
class CcdInit(MSONable, ToJsonFileMixIn):
    name: str
    excited_structure: Structure  # excited excited structure
    ground_structure: Structure  # final ground structure
    excited_charge: int
    ground_charge: int
    excited_energy: float
    excited_energy_correction: float
    ground_energy: float
    ground_energy_correction: float

    @property
    def dQ(self):
        return get_dQ(self.excited_structure, self.ground_structure)

    def __str__(self):
        charge_diff = self.ground_charge - self.excited_charge
        if charge_diff == -1:
            trapped_carrier = "e-"
        elif charge_diff == 1:
            trapped_carrier = "h+"
        else:
            raise ValueError("The charge difference between excited and ground "
                             "states is neither 1 nor -1")

        result = []
        excited_state = f"{self.name}_{self.excited_charge} + {trapped_carrier}"
        result.append(["Excited state:", excited_state,
                       "energy:", self.excited_energy,
                       "correction:", self.excited_energy_correction])

        ground_state = f"{self.name}_{self.ground_charge}"
        result.append(["Ground state:", ground_state,
                       "energy:", self.ground_energy,
                       "correction:", self.ground_energy_correction])

        return tabulate(result, tablefmt="plain")


@dataclass
class ImageStructureInfo(MSONable):
    displace_ratio: float  # 0.0: original, 1.0: counter structure
    energy: float  # include the corrections


@dataclass
class Ccd(MSONable, ToJsonFileMixIn):
    dQ: float
    excited_image_infos: List[ImageStructureInfo]
    ground_image_infos: List[ImageStructureInfo]
    correction: str = None

    @property
    def ground_dQs(self):
        return [self.dQ * i.displace_ratio for i in self.ground_image_infos]

    @property
    def ground_energies(self):
        return [i.energy for i in self.ground_image_infos]

    @property
    def excited_dQs(self):
        return [self.dQ * (1 - i.displace_ratio)
                for i in self.excited_image_infos]

    @property
    def excited_energies(self):
        return [i.energy for i in self.excited_image_infos]


def ccd_plt(ccd: Ccd):
    plt.plot(ccd.ground_dQs, ccd.ground_energies)
    plt.plot(ccd.excited_dQs, ccd.excited_energies)
    return plt

