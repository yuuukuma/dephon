# -*- coding: utf-8 -*-
#  Copyright (c) 2022 Kumagai group.
import os
from argparse import Namespace
from pathlib import Path

from monty.serialization import loadfn
from pydefect.analyzer.defect_energy import DefectEnergyInfo

from dephon.make_config_coord import make_ccd_init


def make_ccd_init_and_dirs(args: Namespace):
    i_calc_results = loadfn(args.initial_dir / "calc_results.json")
    f_calc_results = loadfn(args.final_dir / "calc_results.json")

    i_defect_energy_info = DefectEnergyInfo.from_yaml(
        args.initial_dir / "defect_energy_info.yaml")
    f_defect_energy_info = DefectEnergyInfo.from_yaml(
        args.final_dir / "defect_energy_info.yaml")

    ccd_init = make_ccd_init(i_calc_results, f_calc_results,
                             i_defect_energy_info, f_defect_energy_info,
                             i_to_f_div_ratios=args.i_to_f_div_ratios,
                             f_to_i_div_ratios=args.f_to_i_div_ratios)

    i_charge, f_charge = ccd_init.initial_charge, ccd_init.final_charge
    path = Path(f"cc/{ccd_init.name}_{i_charge}to{f_charge}")
    path.mkdir(parents=True)

    for i, imag_structures in [("initial", ccd_init.i_to_f_image_structures),
                               ("final", ccd_init.f_to_i_image_structures)]:
        os.mkdir(path / i)
        for imag_structure in imag_structures:
            dir_ = path / i / str(imag_structure.displace_ratio)
            os.mkdir(dir_)
            imag_structure.structure.to(filename=str(dir_ / "POSCAR"))

    return ccd_init
