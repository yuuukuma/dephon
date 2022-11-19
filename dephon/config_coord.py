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
    initial_structure: Structure  # initial excited structure
    final_structure: Structure  # final ground structure
    initial_charge: int
    final_charge: int
    initial_energy: DefectEnergy
    final_energy: DefectEnergy
    i_to_f_image_structures: List[ImageStructure]
    f_to_i_image_structures: List[ImageStructure]

    @property
    def dQ(self):
        return get_dQ(self.initial_structure, self.final_structure)

    def __str__(self):
        charge_diff = self.final_charge - self.initial_charge
        if charge_diff == -1:
            trapped_carrier = "e-"
        elif charge_diff == 1:
            trapped_carrier = "h+"
        else:
            raise ValueError("The charge difference between initial and final "
                             "states is neither 1 nor -1")

        initial = f"{self.name}_{self.initial_charge} + {trapped_carrier}"
        final = f"{self.name}_{self.final_charge}"

        return f"""Name: {self.name}
transition: {initial} -> {final}"""


@dataclass
class Ccd(MSONable):
    name: str
    dQs: List[float]
    initial_energies: List[float]
    final_energies: List[float]


