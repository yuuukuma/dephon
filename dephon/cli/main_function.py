# -*- coding: utf-8 -*-
#  Copyright (c) 2022 Kumagai group.
from argparse import Namespace
from glob import glob
from pathlib import Path

import yaml
from monty.serialization import loadfn
from pydefect.analyzer.calc_results import CalcResults
from vise.input_set.prior_info import PriorInfo
from vise.util.logger import get_logger

from dephon.config_coord import ImageStructureInfo, Ccd, ccd_plt, CcdPlotter
from dephon.make_config_coord import make_ccd_init


logger = get_logger(__name__)


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

    ccd_init.to_json_file(str(path / "ccd_init.json"))
    print(ccd_init)

    _make_ccd_dirs(args, ccd_init, path)


def _make_ccd_dirs(args, ccd_init, path: Path):
    gs, es = ccd_init.ground_structure, ccd_init.excited_structure
    e_to_g = es.interpolate(gs, nimages=args.e_to_g_div_ratios)
    g_to_e = gs.interpolate(es, nimages=args.g_to_e_div_ratios)

    for i, ratios, ss in [("excited", args.e_to_g_div_ratios, e_to_g),
                          ("ground", args.g_to_e_div_ratios, g_to_e)]:
        (path / i).mkdir(parents=True, exist_ok=True)
        c = ccd_init.ground_charge if i == "ground" else ccd_init.excited_charge

        for ratio, s in zip(ratios, ss):

            dir_ = path / i / f"disp_{ratio}"
            try:
                dir_.mkdir()
                logger.info(f"Directory {dir_} was created.")
                s.to(filename=str(dir_ / "POSCAR"))
                (dir_ / "prior_info.yaml").write_text(
                    yaml.dump({"charge": c}), None)
            except FileExistsError:
                logger.info(f"Directory {dir_} exists, so skip it.")


def add_ccd_dirs(args: Namespace):
    _make_ccd_dirs(args, args.ccd_init, args.calc_dir)


def make_ccd(args: Namespace):
    dQ = args.ccd_init.dQ

    def _make_image_info(relaxed_structure_energy, dir_name):
        _dQ = 0.0 if dir_name == "ground" else dQ

        result = [ImageStructureInfo(0.0, relaxed_structure_energy, _dQ)]
        for d in glob(f'{dir_name}/disp_*'):
            disp_ratio = float(d.split("_")[-1])
            try:
                cr: CalcResults = loadfn(Path(d) / "calc_results.json")
            except FileNotFoundError:

                continue
            _dQ = dQ * disp_ratio if dir_name == "ground" else dQ * (1. - disp_ratio)
            result.append(ImageStructureInfo(disp_ratio, cr.energy, _dQ))
        result.sort(key=lambda x: x.displace_ratio)
        return result

    g_energy = args.ccd_init.ground_energy
    ground_image_infos = _make_image_info(g_energy, "ground")
    ground_correction = args.ccd_init.ground_energy_correction
    for i in ground_image_infos:
        i.energy += ground_correction

    e_energy = args.ccd_init.excited_energy
    excited_image_infos = _make_image_info(e_energy, "excited")
    excited_correction = args.ccd_init.excited_energy_correction
    for i in excited_image_infos:
        i.energy += excited_correction

    ccd = Ccd(dQ, excited_image_infos, ground_image_infos,
              correction_type="constant FNV")
    ccd.to_json_file()


def plot_ccd(args: Namespace):
    plotter = CcdPlotter(args.ccd, spline_deg=args.spline_deg)
    plotter.construct_plot()
    plotter.plt.savefig(args.fig_name)
    plotter.plt.show()
