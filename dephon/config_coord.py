# -*- coding: utf-8 -*-
#  Copyright (c) 2022 Kumagai group.
from dataclasses import dataclass
from typing import List

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
    excited_energy: DefectEnergy
    ground_energy: DefectEnergy

    @property
    def dQ(self):
        return get_dQ(self.excited_structure, self.ground_structure)

    def __str__(self):
        result = [["Name:", self.name]]
        charge_diff = self.ground_charge - self.excited_charge
        if charge_diff == -1:
            trapped_carrier = "e-"
        elif charge_diff == 1:
            trapped_carrier = "h+"
        else:
            raise ValueError("The charge difference between excited and ground "
                             "states is neither 1 nor -1")

        excited_state = f"{self.name}_{self.excited_charge} + {trapped_carrier}"
        result.append(["Excited state:", excited_state,
                       "energy:", self.excited_energy.formation_energy,
                       "total correction:", self.excited_energy.total_correction,
                       "is shallow:", self.excited_energy.is_shallow])

        ground_state = f"{self.name}_{self.ground_charge}"
        result.append(["Ground state:", ground_state,
                       "energy:", self.ground_energy.formation_energy,
                       "total correction:", self.ground_energy.total_correction,
                       "is shallow:", self.ground_energy.is_shallow])

        return tabulate(result, tablefmt="plain")


@dataclass
class ImageStructureInfo(MSONable):
    displace_ratio: float  # 0.0: original, 1.0: counter structure
    energy: float


@dataclass
class Ccd(MSONable, ToJsonFileMixIn):
    name: str
    dQ: float
    excited_image_structures: List[ImageStructureInfo]
    ground_image_structures: List[ImageStructureInfo]


