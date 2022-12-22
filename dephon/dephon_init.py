# -*- coding: utf-8 -*-
#  Copyright (c) 2022 Kumagai group.
from copy import deepcopy
from dataclasses import dataclass
from pathlib import Path
from typing import List

import numpy as np
from monty.json import MSONable
from pymatgen.analysis.defects.ccd import get_dQ
from pymatgen.core import Structure
from tabulate import tabulate
from vise.util.logger import get_logger
from vise.util.mix_in import ToJsonFileMixIn
from vise.util.structure_symmetrizer import num_sym_op

from dephon.enum import Carrier

logger = get_logger(__name__)


def get_dR(ground: Structure, excited: Structure) -> float:
    """Summation of atomic displacement distance

    Args:
        ground (Structure): Reference structure
        excited (Structure): Target structure

    Constant deviation needs to be removed by e.g., aligning the farthest atom.

    Returns:
        The Summed atomic displacement distance in float
    """
    return np.sqrt(
        np.sum([x.distance(y) ** 2 for x, y in zip(ground, excited)]))


@dataclass
class MinimumPointInfo(MSONable):
    charge: int
    structure: Structure
    # formation energy at Ef=VBM and chemical potentials being standard states.
    energy: float
    energy_correction: float
    carriers: List[Carrier]
    initial_site_symmetry: str
    final_site_symmetry: str
    # absolute dir
    parsed_dir: str

    @property
    def corrected_energy(self):
        return self.energy + self.energy_correction

    @property
    def degeneracy_by_symmetry_reduction(self):
        initial_num_sym_op = num_sym_op[self.initial_site_symmetry]
        final_num_sym_op = num_sym_op[self.final_site_symmetry]
        return initial_num_sym_op / final_num_sym_op

    @property
    def carriers_str(self):
        return "+".join([str(c) for c in self.carriers])

    @property
    def dir_path(self):
        return Path(self.parsed_dir)


@dataclass
class DephonInit(MSONable, ToJsonFileMixIn):
    name: str
    ground_state: MinimumPointInfo
    excited_state: MinimumPointInfo
    vbm: float
    cbm: float
    supercell_vbm: float
    supercell_cbm: float

    def __post_init__(self):
        if abs(self.ground_state.charge - self.excited_state.charge) != 1:
            logger.warning("Charge difference between ground and excited states"
                           " is not 1. You must know what you are doing.")
        self.ground_state_w_pn = deepcopy(self.ground_state)
        self.ground_state_w_pn.carriers = [Carrier.hole, Carrier.electron]
        self.ground_state_w_pn.energy += self.cbm - self.vbm

    @property
    def band_gap(self):
        return self.cbm - self.vbm

    @property
    def delta_EF(self):
        return self.band_gap if self.semiconductor_type == "n" else 0.0

    @property
    def semiconductor_type(self) -> str:
        trap_charge_by_excited_state = \
            - (self.excited_state.charge - self.ground_state.charge)
        if trap_charge_by_excited_state == 1:
            return "p"
        else:
            return "n"

    @property
    def volume(self):
        assert (self.ground_state.structure.volume
                == self.excited_state.structure.volume)
        return self.ground_state.structure.volume

    @property
    def dQ(self):
        return get_dQ(self.excited_state.structure, self.ground_state.structure)

    @property
    def dR(self):
        return get_dR(self.excited_state.structure, self.ground_state.structure)

    @property
    def modal_mass(self):
        return (self.dQ / self.dR) ** 2

    def __str__(self):
        result = [f"name: {self.name}",
                  f"semiconductor type:  {self.semiconductor_type}-type"]
        table = [["vbm", self.vbm, "supercell vbm", self.supercell_vbm],
                 ["cbm", self.cbm, "supercell cbm", self.supercell_cbm],
                 ["dQ (amu^0.5 Å)", self.dQ],
                 ["dR (Å)", self.dR],
                 ["M (amu)", self.modal_mass]]
        result.append(tabulate(table, tablefmt="plain", floatfmt=".3f"))

        result.append("-" * 60)

        headers = ["q", "carrier", "initial symm", "final symm", "energy",
                   "correction", "corrected energy", "ZPL"]
        table = []

        last_energy = None

        for state in [self.ground_state_w_pn,
                      self.excited_state,
                      self.ground_state]:
            table.append(
                [state.charge, state.carriers_str, state.initial_site_symmetry,
                 state.final_site_symmetry, state.energy,
                 state.energy_correction, state.corrected_energy])
            if last_energy:
                table[-1].append(last_energy - state.corrected_energy)
            last_energy = state.corrected_energy

        result.append(
            tabulate(table, tablefmt="plain", headers=headers, floatfmt=".3f",
                     stralign="center"))
        # result.append(f"ZPL: {self.zero_phonon_line:.3f}")
        return "\n".join(result)

""" 
TODO: 
plot

1. add defect position
2. consider how to handle the small difference of origin.
"""

