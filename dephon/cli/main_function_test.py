# -*- coding: utf-8 -*-
#  Copyright (c) 2022 Kumagai group.
import os
import shutil
from argparse import Namespace
from pathlib import Path

from monty.serialization import loadfn
from pydefect.analyzer.calc_results import CalcResults
from pydefect.analyzer.defect_energy import DefectEnergy
from pymatgen.core import Structure
from vise.input_set.prior_info import PriorInfo

from dephon.cli.main_function import make_ccd_init_and_dirs, make_ccd, \
    add_ccd_dirs
from dephon.config_coord import Ccd, CcdInit, ImageStructureInfo


def test_make_ccd_init_and_dirs(test_files, tmpdir):
    print(tmpdir)
    tmpdir.chdir()
    shutil.copytree(test_files / "Na3AgO2" / "Va_O1_0", tmpdir / "Va_O1_0")
    shutil.copytree(test_files / "Na3AgO2" / "Va_O1_1", tmpdir / "Va_O1_1")

    args = Namespace(excited_dir=Path("Va_O1_0"),
                     ground_dir=Path("Va_O1_1"),
                     e_to_g_div_ratios=[0.1, 0.2],
                     g_to_e_div_ratios=[0.3, 0.4])
    make_ccd_init_and_dirs(args)


def test_add_ccd_dirs(tmpdir, mocker, ground_structure, excited_structure,
                      intermediate_structure):
    print(tmpdir)
    tmpdir.chdir()

    Path("test/excited/disp_1.0").mkdir(parents=True)

    ccd_init = mocker.Mock(spec=CcdInit, autospec=True)
    ccd_init.excited_structure = excited_structure
    ccd_init.ground_structure = ground_structure

    ccd_init.excited_charge = -1
    ccd_init.ground_charge = 0

    args = Namespace(ccd_init=ccd_init,
                     e_to_g_div_ratios=[0.5, 1.0],
                     g_to_e_div_ratios=[1.0],
                     calc_dir=Path("test"))
    add_ccd_dirs(args)
    excited_5 = Path("test/excited") / "disp_0.5"
    actual = Structure.from_file(excited_5 / "POSCAR")
    assert actual == intermediate_structure

    actual = PriorInfo.load_yaml(excited_5 / "prior_info.yaml")
    expected = PriorInfo(charge=-1)
    assert actual == expected

    actual = Structure.from_file(Path("test/ground") / "disp_1.0" / "POSCAR")
    assert actual == excited_structure


def test_make_ccd(tmpdir, mocker, ground_structure):
    print(tmpdir)
    tmpdir.chdir()
    dirs = [-0.2, 0.2]
    for i in ["ground", "excited"]:
        for d in dirs:
            dir_ = Path(f"{i}/disp_{d}")
            dir_.mkdir(parents=True)
            calc_results = CalcResults(ground_structure, energy=d,
                                       magnetization=0.0, potentials=[0.0])
            calc_results.to_json_file(str(dir_ / "calc_results.json"))

    ccd_init = mocker.Mock(spec=CcdInit, autospec=True)
    ccd_init.dQ = 10.0
    ccd_init.excited_energy = 10.0
    ccd_init.ground_energy = 20.0
    ccd_init.excited_energy_correction = -1.0
    ccd_init.ground_energy_correction = -2.0

    args = Namespace(ccd_init=ccd_init, fig_name=None)
    make_ccd(args)
    actual: Ccd = loadfn("ccd.json")
    expected = Ccd(dQ=10.0,
                   excited_image_infos=[ImageStructureInfo(-0.2, -1.2),
                                        ImageStructureInfo(0.0, 9.0),
                                        ImageStructureInfo(0.2, -0.8)],
                   ground_image_infos=[ImageStructureInfo(-0.2, -2.2),
                                       ImageStructureInfo(0.0, 18.0),
                                       ImageStructureInfo(0.2, -1.8)],
                   correction="constant FNV")
    assert actual == expected


"""
TODO:
. Check if electronic SCF are converged.
. Check if displace_ratio=1 structure is the same as the counterpart.
  Consider the energy correction
"""