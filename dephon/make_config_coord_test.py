# -*- coding: utf-8 -*-
#  Copyright (c) 2022 Kumagai group.

from pydefect.analyzer.calc_results import CalcResults
from pydefect.analyzer.defect_energy import DefectEnergy, DefectEnergyInfo

from dephon.config_coord import CcdInit, ImageStructure
from dephon.make_config_coord import make_ccd_init


def test_make_cc_init(
        initial_structure, final_structure, mocker, intermediate_structure):

    initial_calc_results = mocker.Mock(spec=CalcResults, autospec=True)
    initial_calc_results.structure = initial_structure
    final_calc_results = mocker.Mock(spec=CalcResults, autospec=True)
    final_calc_results.structure = final_structure

    initial_energy = mocker.Mock(spec=DefectEnergy, autospec=True)
    final_energy = mocker.Mock(spec=DefectEnergy, autospec=True)

    initial_defect_energy_info = DefectEnergyInfo(
        "Va_O1", 1, {"O": -1}, initial_energy)
    final_defect_energy_info = DefectEnergyInfo(
        "Va_O1", 0, {"O": -1}, final_energy)

    actual = make_ccd_init(initial_calc_results,
                           final_calc_results,
                           initial_defect_energy_info,
                           final_defect_energy_info,
                           [0.5], [0.5])
    expected = CcdInit("Va_O1",
                       initial_structure=initial_structure,
                       final_structure=final_structure,
                       initial_charge=1,
                       final_charge=0,
                       initial_energy=initial_energy,
                       final_energy=final_energy,
                       i_to_f_image_structures=
                       [ImageStructure(intermediate_structure, 0.5)],
                       f_to_i_image_structures=
                       [ImageStructure(intermediate_structure, 0.5)])
    assert actual == expected
