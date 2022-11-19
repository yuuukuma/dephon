# -*- coding: utf-8 -*-
#  Copyright (c) 2022 Kumagai group.
from typing import List

from pydefect.analyzer.calc_results import CalcResults
from pydefect.analyzer.defect_energy import DefectEnergyInfo

from dephon.config_coord import CcdInit, ImageStructure


def make_ccd_init(initial_calc_results: CalcResults,
                  final_calc_results: CalcResults,
                  initial_defect_energy_info: DefectEnergyInfo,
                  final_defect_energy_info: DefectEnergyInfo,
                  i_to_f_div_ratios: List[float],
                  f_to_i_div_ratios: List[float]) -> CcdInit:
    initial_charge = initial_defect_energy_info.charge
    final_charge = final_defect_energy_info.charge

    assert abs(initial_charge - final_charge) == 1
    assert initial_defect_energy_info.name == final_defect_energy_info.name
    assert (initial_defect_energy_info.atom_io
            == final_defect_energy_info.atom_io)

    i_structure = initial_calc_results.structure
    f_structure = final_calc_results.structure

    i_to_f = i_structure.interpolate(f_structure, nimages=i_to_f_div_ratios)
    f_to_i = f_structure.interpolate(i_structure, nimages=f_to_i_div_ratios)

    i_to_f_s = [ImageStructure(s, d) for s, d in zip(i_to_f, i_to_f_div_ratios)]
    f_to_i_s = [ImageStructure(s, d) for s, d in zip(f_to_i, f_to_i_div_ratios)]

    return CcdInit(name=initial_defect_energy_info.name,
                   initial_structure=i_structure,
                   final_structure=f_structure,
                   initial_charge=initial_charge,
                   final_charge=final_charge,
                   initial_energy=initial_defect_energy_info.defect_energy,
                   final_energy=final_defect_energy_info.defect_energy,
                   i_to_f_image_structures=i_to_f_s,
                   f_to_i_image_structures=f_to_i_s)
