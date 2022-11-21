# -*- coding: utf-8 -*-
#  Copyright (c) 2022 Kumagai group.
from typing import List

from pydefect.analyzer.calc_results import CalcResults
from pydefect.analyzer.defect_energy import DefectEnergyInfo

from dephon.config_coord import CcdInit, ImageStructureInfo


def make_ccd_init(excited_calc_results: CalcResults,
                  ground_calc_results: CalcResults,
                  excited_defect_energy_info: DefectEnergyInfo,
                  ground_defect_energy_info: DefectEnergyInfo) -> CcdInit:
    excited_charge = excited_defect_energy_info.charge
    ground_charge = ground_defect_energy_info.charge

    assert abs(excited_charge - ground_charge) == 1
    assert excited_defect_energy_info.name == ground_defect_energy_info.name
    assert (excited_defect_energy_info.atom_io
            == ground_defect_energy_info.atom_io)

    e_structure = excited_calc_results.structure
    g_structure = ground_calc_results.structure

    return CcdInit(name=excited_defect_energy_info.name,
                   excited_structure=e_structure,
                   ground_structure=g_structure,
                   excited_charge=excited_charge,
                   ground_charge=ground_charge,
                   excited_energy=excited_defect_energy_info.defect_energy,
                   ground_energy=ground_defect_energy_info.defect_energy)
