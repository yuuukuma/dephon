# -*- coding: utf-8 -*-
#  Copyright (c) 2022 Kumagai group.
from dataclasses import dataclass
from typing import List

from monty.json import MSONable
from pydefect.analyzer.defect_energy import DefectEnergy
from pymatgen.analysis.defects.ccd import get_dQ
from pymatgen.core import Structure
from vise.util.mix_in import ToJsonFileMixIn


@dataclass
class ImageStructure(MSONable):
    structure: Structure
    displace_ratio: float  # 0.0: original, 1.0: counter structure


@dataclass
class CcdInit(MSONable, ToJsonFileMixIn):
    name: str
    excited_structure: Structure  # excited excited structure
    ground_structure: Structure  # final ground structure
    excited_charge: int
    ground_charge: int
    excited_energy: DefectEnergy
    ground_energy: DefectEnergy
    e_to_g_image_structures: List[ImageStructure]
    g_to_e_image_structures: List[ImageStructure]

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

        excited = f"{self.name}_{self.excited_charge} + {trapped_carrier}"
        ground = f"{self.name}_{self.ground_charge}"

        return f"""Name: {self.name}
transition: {excited} -> {ground},
"""


@dataclass
class Ccd(MSONable):
    name: str
    dQs: List[float]
    excited_energies: List[float]
    ground_energies: List[float]


