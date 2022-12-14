# -*- coding: utf-8 -*-
#  Copyright (c) 2022 Kumagai group.
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
    """ Information at lowest energy point at given charge

    Attributes:
        charge (int): The charge state.
        structure (Structure): The atomic configuration.
        energy (float): Formation energy at Ef=VBM and chemical potentials
            being standard states.
        correction_energy (float): Correction energy estimated e.g. eFNV.
        initial_site_symmetry (str): Site symmetry before relaxing a defect.
        final_site_symmetry (str): Site symmetry after relaxing a defect.
        parse_dir (str): Directory where the calculation results of this
            minimum point are stored. This should be an absolute path.
    """
    charge: int
    structure: Structure
    energy: float
    correction_energy: float
    initial_site_symmetry: str
    final_site_symmetry: str
    # absolute dir
    parsed_dir: str

    @property
    def corrected_energy(self):
        return self.energy + self.correction_energy

    @property
    def degeneracy_by_symmetry_reduction(self):
        initial_num_sym_op = num_sym_op[self.initial_site_symmetry]
        final_num_sym_op = num_sym_op[self.final_site_symmetry]
        return initial_num_sym_op / final_num_sym_op

    @property
    def dir_path(self):
        return Path(self.parsed_dir)


@dataclass
class DephonInit(MSONable, ToJsonFileMixIn):
    """ Information related to configuration coordination diagram.

    Attributes:
        defect_name (str): Name of a defect, e.g., Va_O1
        states (List[MinimumPointInfo]): List of two minimum points. Usually,
            the charge state difference should be 1.
        vbm (float): valence band maximum in the unitcell calculation.
        cbm (float): conduction band minimum in the unitcell calculation.
        superell_vbm (float): vbm in the perfect supercell calculation.
        superell_cbm (float): cbm in the perfect supercell calculation.
    """
    defect_name: str
    states: List[MinimumPointInfo]
    vbm: float
    cbm: float
    supercell_vbm: float
    supercell_cbm: float

    def __post_init__(self):
        assert len(self.states) == 2

    @property
    def band_gap(self):
        return self.cbm - self.vbm

    def min_point_info_from_charge(self, charge: int):
        for state in self.states:
            if state.charge == charge:
                return state
        raise ValueError(f"Charge {charge} does not exist.")

    @property
    def volume(self):
        assert (self.states[0].structure.volume
                == self.states[1].structure.volume)
        return self.states[0].structure.volume

    @property
    def dQ(self):
        return abs(get_dQ(self.states[0].structure, self.states[1].structure))

    @property
    def dR(self):
        return abs(get_dR(self.states[0].structure, self.states[1].structure))

    @property
    def modal_mass(self):
        return (self.dQ / self.dR) ** 2

    def __str__(self):
        result = [f"name: {self.defect_name}"]
        table = [["vbm", self.vbm, "supercell vbm", self.supercell_vbm],
                 ["cbm", self.cbm, "supercell cbm", self.supercell_cbm],
                 ["dQ (amu^0.5 ??)", self.dQ],
                 ["dR (??)", self.dR],
                 ["M (amu)", self.modal_mass]]
        result.append(tabulate(table, tablefmt="plain", floatfmt=".3f"))

        result.append("-" * 60)

        headers = ["q", "initial symm", "final symm", "energy",
                   "correction", "corrected energy", "ZPL"]
        table = []

        last_energy = None

        for state in self.states:
            table.append(
                [state.charge, state.initial_site_symmetry,
                 state.final_site_symmetry, state.energy,
                 state.correction_energy, state.corrected_energy])
            if last_energy:
                table[-1].append(last_energy - state.corrected_energy)
            last_energy = state.corrected_energy

        result.append(
            tabulate(table, tablefmt="plain", headers=headers, floatfmt=".3f",
                     stralign="center"))
        return "\n".join(result)

""" 
TODO: 
plot

1. add defect position
2. consider how to handle the small difference of origin.
"""

