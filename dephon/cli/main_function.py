# -*- coding: utf-8 -*-
#  Copyright (c) 2022 Kumagai group.
import os
from argparse import Namespace
from glob import glob
from pathlib import Path

from monty.serialization import loadfn
from pydefect.analyzer.calc_results import CalcResults
from pydefect.analyzer.defect_energy import DefectEnergyInfo

from dephon.config_coord import ImageStructureInfo, Ccd, ccd_plt
from dephon.make_config_coord import make_ccd_init


def make_ccd_init_and_dirs(args: Namespace):
    name_e = "_".join(args.excited_dir.name.split("_")[:-1])
    name_g = "_".join(args.ground_dir.name.split("_")[:-1])
    assert name_e == name_g

    e_calc_results = loadfn(args.excited_dir / "calc_results.json")
    g_calc_results = loadfn(args.ground_dir / "calc_results.json")

    e_correction = loadfn(args.excited_dir / "correction.json")
    g_correction = loadfn(args.ground_dir / "correction.json")

    ccd_init = make_ccd_init(name_e, e_calc_results, g_calc_results,
                             e_correction, g_correction)

    e_charge, g_charge = ccd_init.excited_charge, ccd_init.ground_charge
    path = Path(f"cc/{ccd_init.name}_{e_charge}to{g_charge}")
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

    ccd_init.to_json_file(str(path / "ccd_init.json"))
    print(ccd_init)


def _make_image_info(relaxed_structure_energy, dir_name):
    result = [ImageStructureInfo(0.0, energy=relaxed_structure_energy)]
    for d in glob(f'{dir_name}/disp_*'):
        disp_ratio = float(d.split("_")[-1])
        cr: CalcResults = loadfn(Path(d) / "calc_results.json")
        result.append(ImageStructureInfo(disp_ratio, cr.energy))
    result.sort(key=lambda x: x.displace_ratio)
    return result


def make_ccd(args: Namespace):
    g_energy = args.ccd_init.ground_energy.energy(False)
    ground_image_infos = _make_image_info(g_energy, "ground")
    ground_correction = args.ccd_init.ground_energy.total_correction
    for i in ground_image_infos:
        i.energy += ground_correction

    e_energy = args.ccd_init.excited_energy.energy(False)
    excited_image_infos = _make_image_info(e_energy, "excited")
    excited_correction = args.ccd_init.excited_energy.total_correction
    for i in excited_image_infos:
        i.energy += excited_correction

    ccd = Ccd(args.ccd_init.dQ, excited_image_infos, ground_image_infos,
              correction="constant FNV")
    ccd.to_json_file()
    plt = ccd_plt(ccd)
    plt.show()
    plt.savefig("CCD.pdf")