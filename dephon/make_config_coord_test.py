# -*- coding: utf-8 -*-
#  Copyright (c) 2022 Kumagai group.
from pydefect.analyzer.band_edge_states import PerfectBandEdgeState, EdgeInfo, \
    OrbitalInfo
from pydefect.analyzer.calc_results import CalcResults
from pydefect.analyzer.unitcell import Unitcell
from pydefect.corrections.efnv_correction import ExtendedFnvCorrection

from dephon.config_coord import CcdInit
from dephon.make_config_coord import make_ccd_init


def test_make_cc_init(
        excited_structure, ground_structure, mocker):

    excited_calc_results = mocker.Mock(spec=CalcResults, autospec=True)
    excited_calc_results.structure = excited_structure
    excited_calc_results.energy = 10.0

    ground_calc_results = mocker.Mock(spec=CalcResults, autospec=True)
    ground_calc_results.structure = ground_structure
    ground_calc_results.energy = 1.0

    e_correction = mocker.Mock(spec=ExtendedFnvCorrection, autospec=True)
    e_correction.correction_energy = 100.0
    e_correction.charge = 1

    g_correction = mocker.Mock(spec=ExtendedFnvCorrection, autospec=True)
    g_correction.correction_energy = 1000.0
    g_correction.charge = 0

    unitcell = mocker.Mock(spec=Unitcell, autospec=True)
    unitcell.vbm = 0.0
    unitcell.cbm = 10.0

    p_band_edge_state = mocker.Mock(spec=PerfectBandEdgeState, autospec=True)
    p_band_edge_state.vbm_info = EdgeInfo(1, (0.0, 0.0, 0.0),
                                          OrbitalInfo(2.0, {}, 0.0))
    p_band_edge_state.cbm_info = EdgeInfo(1, (0.0, 0.0, 0.0),
                                          OrbitalInfo(8.0, {}, 0.0))

    actual = make_ccd_init("Va_O1",
                           excited_calc_results,
                           ground_calc_results,
                           e_correction,
                           g_correction,
                           unitcell,
                           p_band_edge_state)
    expected = CcdInit("Va_O1",
                       excited_structure=excited_structure,
                       ground_structure=ground_structure,
                       excited_charge=1,
                       ground_charge=0,
                       excited_energy=10.0,
                       excited_energy_correction=100.0,
                       ground_energy=1.0,
                       ground_energy_correction=1000.0,
                       vbm=0.0, cbm=10.0,
                       supercell_vbm=2.0, supercell_cbm=8.0)
    assert actual == expected
