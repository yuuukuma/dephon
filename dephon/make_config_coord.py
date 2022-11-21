# -*- coding: utf-8 -*-
#  Copyright (c) 2022 Kumagai group.

from pydefect.analyzer.calc_results import CalcResults
from pydefect.corrections.efnv_correction import ExtendedFnvCorrection

from dephon.config_coord import CcdInit


def make_ccd_init(name: str,
                  excited_calc_results: CalcResults,
                  ground_calc_results: CalcResults,
                  excited_efnv_corr: ExtendedFnvCorrection,
                  ground_efnv_corr: ExtendedFnvCorrection) -> CcdInit:

    assert abs(excited_efnv_corr.charge - ground_efnv_corr.charge) == 1
    return CcdInit(
        name=name,
        excited_structure=excited_calc_results.structure,
        ground_structure=ground_calc_results.structure,
        excited_charge=excited_efnv_corr.charge,
        ground_charge=ground_efnv_corr.charge,
        excited_energy=excited_calc_results.energy,
        excited_energy_correction=excited_efnv_corr.correction_energy,
        ground_energy_correction=ground_efnv_corr.correction_energy,
        ground_energy=ground_calc_results.energy)
