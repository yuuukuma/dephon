# -*- coding: utf-8 -*-
#  Copyright (c) 2022 Kumagai group.
import os
from argparse import Namespace
from pathlib import Path

from monty.serialization import loadfn
from pydefect.analyzer.defect_energy import DefectEnergyInfo

from dephon.make_config_coord import make_ccd_init


def make_ccd_init_and_dirs(args: Namespace):
    i_calc_results = loadfn(args.excited_dir / "calc_results.json")
    f_calc_results = loadfn(args.ground_dir / "calc_results.json")

    i_defect_energy_info = DefectEnergyInfo.from_yaml(
        args.excited_dir / "defect_energy_info.yaml")
    f_defect_energy_info = DefectEnergyInfo.from_yaml(
        args.ground_dir / "defect_energy_info.yaml")

    ccd_init = make_ccd_init(i_calc_results, f_calc_results,
                             i_defect_energy_info, f_defect_energy_info)

    i_charge, f_charge = ccd_init.excited_charge, ccd_init.ground_charge
    path = Path(f"cc/{ccd_init.name}_{i_charge}to{f_charge}")
    path.mkdir(parents=True)

    gs, es = ccd_init.ground_structure, ccd_init.excited_structure

    e_to_g = es.interpolate(gs, nimages=args.g_to_e_div_ratios)
    g_to_e = gs.interpolate(es, nimages=args.e_to_g_div_ratios)

    for i, ratios, ss in [("excited", args.e_to_g_div_ratios, e_to_g),
                          ("ground", args.e_to_g_div_ratios, g_to_e)]:
        os.mkdir(path / i)
        for ratio, s in zip(ratios, ss):
            dir_ = path / i / f"disp_{ratio}"
            os.mkdir(dir_)
            s.to(filename=str(dir_ / "POSCAR"))

        if Path(path / i / "disp_0.0").is_dir() is False:
            if i == "excited":
                os.symlink(f"../../../{args.excited_dir}", path / i / "disp_0.0")
            if i == "ground":
                os.symlink(f"../../../{args.ground_dir}", path / i / "disp_0.0")
    ccd_init.to_json_file(str(path / "ccd_init.json"))
    print(ccd_init)
