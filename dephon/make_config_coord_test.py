# -*- coding: utf-8 -*-
#  Copyright (c) 2022 Kumagai group.

from pydefect.analyzer.calc_results import CalcResults
from pydefect.analyzer.defect_energy import DefectEnergy, DefectEnergyInfo

from dephon.config_coord import CcdInit
from dephon.make_config_coord import make_ccd_init


def test_make_cc_init(
        excited_structure, ground_structure, mocker):

    excited_calc_results = mocker.Mock(spec=CalcResults, autospec=True)
    excited_calc_results.structure = excited_structure
    ground_calc_results = mocker.Mock(spec=CalcResults, autospec=True)
    ground_calc_results.structure = ground_structure

    excited_energy = mocker.Mock(spec=DefectEnergy, autospec=True)
    ground_energy = mocker.Mock(spec=DefectEnergy, autospec=True)

    excited_defect_energy_info = DefectEnergyInfo(
        "Va_O1", 1, {"O": -1}, excited_energy)
    ground_defect_energy_info = DefectEnergyInfo(
        "Va_O1", 0, {"O": -1}, ground_energy)

    actual = make_ccd_init(excited_calc_results,
                           ground_calc_results,
                           excited_defect_energy_info,
                           ground_defect_energy_info)
    expected = CcdInit("Va_O1",
                       excited_structure=excited_structure,
                       ground_structure=ground_structure,
                       excited_charge=1,
                       ground_charge=0,
                       excited_energy=excited_energy,
                       ground_energy=ground_energy)
    assert actual == expected
