# -*- coding: utf-8 -*-
#  Copyright (c) 2022 Kumagai group.
from typing import List

from pydefect.analyzer.calc_results import CalcResults
from pydefect.analyzer.defect_energy import DefectEnergyInfo

from dephon.config_coord import CcdInit, ImageStructureInfo


def make_ccd_init(excited_calc_results: CalcResults,
                  ground_calc_results: CalcResults,
                  excited_defect_energy_info: DefectEnergyInfo,
                  ground_defect_energy_info: DefectEnergyInfo,
                  e_to_g_div_ratios: List[float],
                  g_to_e_div_ratios: List[float]) -> CcdInit:
    excited_charge = excited_defect_energy_info.charge
    ground_charge = ground_defect_energy_info.charge

    assert abs(excited_charge - ground_charge) == 1
    assert excited_defect_energy_info.name == ground_defect_energy_info.name
    assert (excited_defect_energy_info.atom_io
            == ground_defect_energy_info.atom_io)

    e_structure = excited_calc_results.structure
    g_structure = ground_calc_results.structure

    e_to_g = e_structure.interpolate(g_structure, nimages=e_to_g_div_ratios)
    g_to_e = g_structure.interpolate(e_structure, nimages=g_to_e_div_ratios)

    e_to_g_s = [ImageStructureInfo(s, d) for s, d in zip(e_to_g, e_to_g_div_ratios)]
    g_to_e_s = [ImageStructureInfo(s, d) for s, d in zip(g_to_e, g_to_e_div_ratios)]

    return CcdInit(name=excited_defect_energy_info.name,
                   excited_structure=e_structure,
                   ground_structure=g_structure,
                   excited_charge=excited_charge,
                   ground_charge=ground_charge,
                   excited_energy=excited_defect_energy_info.defect_energy,
                   ground_energy=ground_defect_energy_info.defect_energy,
                   e_to_g_image_structures=e_to_g_s,
                   g_to_e_image_structures=g_to_e_s)
