# -*- coding: utf-8 -*-
#  Copyright (c) 2022 Kumagai group.
from collections import defaultdict
from dataclasses import dataclass
from pathlib import Path
from typing import List

import numpy as np
from monty.json import MSONable
from pydefect.analyzer.band_edge_states import LocalizedOrbital
from pymatgen.analysis.defects.ccd import get_dQ
from pymatgen.core import Structure
from pymatgen.electronic_structure.core import Spin
from tabulate import tabulate
from vise.util.logger import get_logger
from vise.util.mix_in import ToJsonFileMixIn
from vise.util.structure_symmetrizer import num_sym_op

from dephon.enum import Carrier
from dephon.util import spin_to_idx

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
class NearEdgeState(MSONable):
    band_index: int  # begin from 1.
    kpt_coord: List[float]
    kpt_weight: float
    kpt_index: int  # begin from 1.
    eigenvalue: float
    occupation: float

    def __str__(self):
        k_coord = " ".join([f"{x:.2f}" for x in self.kpt_coord])
        k = [f"index : {self.kpt_index}",
             f"coord: {k_coord}",
             f"weight: {self.kpt_weight}"]
        x = [f"band index: {self.band_index}",
             f"kpt info: ({', '.join(k)})",
             f"eigenvalue: {self.eigenvalue:.2f}",
             f"occupation: {self.occupation:.2f}"]
        return ", ".join(x)


@dataclass
class MinimumPointInfo(MSONable):
    """ Information at lowest energy point at given charge

    Attributes:
        charge: The charge state.
        structure: The atomic configuration.
        energy: Formation energy at Ef=VBM and chemical potentials
            being standard states.
        correction_energy: Correction energy estimated e.g. eFNV.
        magnetization: Magnetization
        localized_orbitals: List of localized orbitals at each spin channel.
            [Spin up orbitals, Spin down orbitals]
        initial_site_symmetry (str): Site symmetry before relaxing _default_single_ccd_for_e_p_coupling defect.
        final_site_symmetry (str): Site symmetry after relaxing _default_single_ccd_for_e_p_coupling defect.
        parse_dir (str): Directory where the calculation results of this
            minimum point are stored. This should be an absolute path.
    """
    charge: int
    structure: Structure
    energy: float
    correction_energy: float
    magnetization: float
    localized_orbitals: List[List[LocalizedOrbital]]
    initial_site_symmetry: str
    final_site_symmetry: str
    # absolute dir
    parsed_dir: str
    valence_bands: List[List[NearEdgeState]]  # [spin][bands]
    conduction_bands: List[List[NearEdgeState]]  # [spin][bands]

    # @property
    # def __post_init__(self):
    #     if self.valence_bands is None:
    #         self.valence_bands = [[]]
    #     if self.conduction_bands is None:
    #         self.conduction_bands = [[]]

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

    def near_edge_states(self,
                         captured_carrier: Carrier,
                         spin: int) -> List[NearEdgeState]:
        bands = self.conduction_bands \
            if captured_carrier is Carrier.e else self.valence_bands
        idx = 0 if len(bands) == 1 else spin
        return bands[idx]

    @property
    def relevant_band_indices(self):
        result = defaultdict(list)
        for spin in [Spin.up, Spin.down]:
            s_idx = spin_to_idx(spin, count_from_1=True)
            kpt_indices = []
            if self.valence_bands:
                for state in self.valence_bands[spin_to_idx(spin)]:
                    result[(s_idx, state.kpt_index)].append(state.band_index)
                    kpt_indices.append(state.kpt_index)
            for state in self.conduction_bands[spin_to_idx(spin)]:
                result[(s_idx, state.kpt_index)].append(state.band_index)
                kpt_indices.append(state.kpt_index)
            for lo in self.localized_orbitals[spin_to_idx(spin)]:
                for kpt_index in kpt_indices:
                    if lo.band_idx not in result[(s_idx, kpt_index)]:
                        result[(s_idx, kpt_index)].append(lo.band_idx)
        return dict(result)


def make_near_edge_states_from_band_edge_orbital_infos():
    pass


@dataclass
class DephonInit(MSONable, ToJsonFileMixIn):
    """ Information related to configuration coordination diagram.

    Attributes:
        defect_name (str): Name of _default_single_ccd_for_e_p_coupling defect, e.g., Va_O1
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
    ave_electron_mass: float
    ave_hole_mass: float
    ave_static_diele_const: float

    def __post_init__(self):
        assert len(self.states) == 2

    @property
    def band_gap(self):
        return self.cbm - self.vbm

    def min_info_from_charge(self, charge: int):
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
                 ["dQ (amu^0.5 Å)", self.dQ],
                 ["dR (Å)", self.dR],
                 ["M (amu)", self.modal_mass],
                 ["electron mass (m0)", self.ave_electron_mass],
                 ["hole mass (m0)", self.ave_hole_mass],
                 ["static diele", self.ave_static_diele_const]]
        result.append(tabulate(table, tablefmt="plain", floatfmt=".3f"))

        result.append("-" * 60)

        headers = ["q", "ini symm", "final symm", "energy",
                   "correction", "corrected energy", "magnetization",
                   "localized state idx", "ZPL"]
        table = []

        last_energy = None

        for state in self.states:

            localized_state_idxs = []
            for s, spin in zip(state.localized_orbitals, ["up", "down"]):
                for ss in s:
                    localized_state_idxs.append(f"{spin}-{ss.band_idx}")
            table.append(
                [state.charge, state.initial_site_symmetry,
                 state.final_site_symmetry, state.energy,
                 state.correction_energy, state.corrected_energy,
                 state.magnetization,
                 ", ".join(localized_state_idxs)])
            if last_energy:
                table[-1].append(last_energy - state.corrected_energy)
            last_energy = state.corrected_energy

        result.append(
            tabulate(table, tablefmt="plain", headers=headers, floatfmt=".3f",
                     stralign="center"))

        for min_info in self.states:
            result.append(f"- q={min_info.charge}")
            for vb, spin in zip(min_info.valence_bands, ["up", "down"]):
                result.append(f"-- valence bands, spin-{spin}")
                result.extend([f"{x}" for x in vb])
            result.append(f"")
            for cb, spin in zip(min_info.conduction_bands, ["up", "down"]):
                result.append(f"-- conduction bands, spin-{spin}")
                result.extend([f"{x}" for x in cb])
            result.append(f"")

        return "\n".join(result)

""" 
TODO: 
plot

1. add defect position
2. consider how to handle the small difference of origin.
"""

