# -*- coding: utf-8 -*-
#  Copyright (c) 2022 Kumagai group.
import os
import shutil
from argparse import Namespace
from pathlib import Path

import pytest
from monty.serialization import loadfn
from pydefect.analyzer.band_edge_states import PerfectBandEdgeState, EdgeInfo, \
    OrbitalInfo
from pydefect.analyzer.calc_results import CalcResults
from pydefect.analyzer.defect_energy import DefectEnergy
from pydefect.analyzer.unitcell import Unitcell
from pydefect.corrections.efnv_correction import ExtendedFnvCorrection
from pymatgen.core import Structure
from vise.input_set.prior_info import PriorInfo

from dephon.cli.main_function import make_ccd_init, make_ccd, plot_ccd, \
    make_ccd_dirs
from dephon.config_coord import Ccd, CcdInit, ImageStructureInfo


def test_add_ccd_dirs(test_files, tmpdir, mocker):
    tmpdir.chdir()
    unitcell = mocker.Mock(spec=Unitcell, autospec=True)
    unitcell.vbm = 0.0
    unitcell.cbm = 10.0

    p_band_edge_state = mocker.Mock(spec=PerfectBandEdgeState, autospec=True)
    p_band_edge_state.vbm_info = EdgeInfo(1, (0.0, 0.0, 0.0),
                                      OrbitalInfo(2.0, {}, 0.0))
    p_band_edge_state.cbm_info = EdgeInfo(1, (0.0, 0.0, 0.0),
                                      OrbitalInfo(8.0, {}, 0.0))

    _dir = test_files / "Na3AgO2"

    args = Namespace(excited_dir=_dir / "Va_O1_0",
                     ground_dir=_dir / "Va_O1_1",
                     unitcell=unitcell,
                     perfect_band_edge_state=p_band_edge_state)
    make_ccd_init(args)
    actual = loadfn("cc/Va_O1_0to1/ccd_init.json").__str__()
    expected = """dQ              2.24
dR              0.32
M               48.37
Excited state:  Va_O1_0 + h+  energy:  -458.348  correction:  0
Ground state:   Va_O1_1       energy:  -461.49   correction:  0.238665"""
    assert actual == expected


def test_make_ccd_dirs(tmpdir, mocker, ground_structure, excited_structure, 
                       intermediate_structure):
    tmpdir.chdir()
    ccd_init = mocker.Mock(spec=CcdInit, autospec=True)
    ccd_init.dQ = 10.0
    ccd_init.ground_structure = ground_structure
    ccd_init.excited_structure = excited_structure
    ccd_init.ground_correction = 100.0
    ccd_init.excited_correction = 200.0

    ccd_init.excited_charge = -1
    ccd_init.ground_charge = 0

    Path("test").mkdir()

    args = Namespace(ccd_init=ccd_init,
                     e_to_g_div_ratios=[0.5, 1.0],
                     g_to_e_div_ratios=[1.0],
                     calc_dir=Path("test"))
    make_ccd_dirs(args)

    actual = Structure.from_file("excited/disp_0.5/POSCAR")
    assert actual == intermediate_structure

    actual = PriorInfo.load_yaml("excited/disp_0.5/prior_info.yaml")
    expected = PriorInfo(charge=-1)
    assert actual == expected

    actual = loadfn("excited/disp_0.5/image_structure_info.json")
    expected = ImageStructureInfo(dQ=5.0, correction=200.0, correction_type="eFNV")
    assert actual == expected

    actual = Structure.from_file("ground/disp_1.0/POSCAR")
    assert actual == excited_structure


# @pytest.fixture
# def ccd():
#     return Ccd(name="test",
#                image_infos=[ImageStructureInfo(-0.2, -1.2, 12.0),
#                                     ImageStructureInfo(0.0, 9.0, 10.0),
#                                     ImageStructureInfo(0.2, -0.8, 8.0)],
#                ground_image_infos=[ImageStructureInfo(-0.2, -2.2, -2.0),
#                                    ImageStructureInfo(0.0, 18.0, 0.0),
#                                    ImageStructureInfo(0.2, -1.8, 2.0)])


def test_make_ccd(tmpdir, mocker, ground_structure):
    print(tmpdir)
    tmpdir.chdir()

    for i in ["ground", "excited"]:
        for disp in [-0.2, 0.2]:
            disp_dir = Path(f"{i}/disp_{disp}")
            disp_dir.mkdir(parents=True)
            image_structure_info = ImageStructureInfo(dQ=10.0 + disp,
                                                      correction=0.1,
                                                      correction_type="eFNV")
            image_structure_info.to_json_file(disp_dir / "image_structure_info.json")
            calc_results = CalcResults(ground_structure, energy=disp,
                                       magnetization=0.0, potentials=[0.0])
            calc_results.to_json_file(str(disp_dir / "calc_results.json"))

    ccd_init = mocker.Mock(spec=CcdInit, autospec=True)
    ccd_init.name = "test"
    ccd_init.vbm = 1.0
    ccd_init.cbm = 2.0
    ccd_init.ground_charge = 1
    ccd_init.excited_charge = 0

    args = Namespace(ccd_init=ccd_init,
                     ground_dirs=[Path("ground/disp_-0.2"), Path("ground/disp_0.2")],
                     excited_dirs=[Path("excited/disp_-0.2"), Path("excited/disp_0.2")])
    print(Path.cwd())
    make_ccd(args)
    actual: Ccd = loadfn("ccd.json")

    ground_info = [ImageStructureInfo(9.8, -0.2 + 1.0, 0.1, "eFNV"),
                   ImageStructureInfo(10.2, 0.2 + 1.0, 0.1, "eFNV")]
    excited_info = [ImageStructureInfo(9.8, -0.2, 0.1, "eFNV"),
                    ImageStructureInfo(10.2, 0.2, 0.1, "eFNV")]
    ground_info_eg = [ImageStructureInfo(9.8, -0.2 + 2.0, 0.1, "eFNV"),
                      ImageStructureInfo(10.2, 0.2 + 2.0, 0.1, "eFNV")]
    expected = Ccd(name="test",
                   image_infos={"ground": ground_info,
                                "excited + p": excited_info,
                                "ground + p + n": ground_info_eg})
    assert actual == expected


def test_plot_ccd(ccd, tmpdir):
    args = Namespace(ccd=ccd, spline_deg=1, fig_name=tmpdir / "ccd.pdf")
    plot_ccd(args)


# def test_plot_eigenvalues(tmpdir, mocker, ground_structure):
#     args = Namespace(dirs=["a"],
#                      ccd=)




"""
TODO:
. Check if electronic SCF are converged.
. Check if displace_ratio=1 structure is the same as the counterpart.
  Consider the energy correction
"""