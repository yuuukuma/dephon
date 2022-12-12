# -*- coding: utf-8 -*-
#  Copyright (c) 2022 Kumagai group.
import os
from argparse import Namespace
from copy import deepcopy
from glob import glob
from pathlib import Path

import yaml
from monty.serialization import loadfn
from pydefect.analyzer.band_edge_states import PerfectBandEdgeState
from pydefect.analyzer.calc_results import CalcResults
from pydefect.analyzer.defect_energy import DefectEnergyInfo
from pydefect.analyzer.unitcell import Unitcell
from pydefect.cli.main_functions import get_calc_results
from pydefect.cli.main_tools import parse_dirs
from pydefect.corrections.efnv_correction import ExtendedFnvCorrection
from vise.util.logger import get_logger

from dephon.config_coord import ImageStructureInfo, Ccd, CcdPlotter, CcdInit
from dephon.plot_eigenvalues import EigenvaluePlotter

logger = get_logger(__name__)


def _parse_dir(_dir):
    energy_info = DefectEnergyInfo.from_yaml(_dir / "defect_energy_info.yaml")
    calc_results = loadfn(_dir / "calc_results.json")
    correction = loadfn(_dir / "correction.json")
    return energy_info, calc_results, correction


def _check_charge_difference(e_charge, g_charge):
    if abs(g_charge - e_charge) != 1:
        logger.critical(
            f"Ground and excited charge states are {g_charge} and {e_charge}. "
            "The difference is not 1. You should understand what you're doing.")


def _check_ground_exited_names(e_name, g_name):
    if g_name != e_name:
        logger.warning("The names of ground and excited states are "
                       f"{g_name} and {e_name}. Here, {g_name} is used.")


def make_ccd_init(args: Namespace):
    g_energy_info, g_calc_results, g_correction = _parse_dir(args.ground_dir)
    e_energy_info, e_calc_results, e_correction = _parse_dir(args.excited_dir)

    _check_charge_difference(e_energy_info.charge, g_energy_info.charge)
    _check_ground_exited_names(e_energy_info.name, g_energy_info.name)

    ccd_init = CcdInit(
        name=g_energy_info.name,
        excited_structure=e_calc_results.structure,
        ground_structure=g_calc_results.structure,
        excited_charge=e_correction.charge,
        ground_charge=g_correction.charge,
        excited_energy=e_calc_results.energy,
        ground_energy=g_calc_results.energy,
        excited_energy_correction=e_correction.correction_energy,
        ground_energy_correction=g_correction.correction_energy,
        vbm=args.unitcell.vbm,
        cbm=args.unitcell.cbm,
        supercell_vbm=args.perfect_band_edge_state.vbm_info.energy,
        supercell_cbm=args.perfect_band_edge_state.cbm_info.energy)

    transfer_name = f"{ccd_init.excited_charge}to{ccd_init.ground_charge}"
    path = Path(f"cc/{ccd_init.name}_{transfer_name}")
    path.mkdir(parents=True)

    ccd_init.to_json_file(str(path / "ccd_init.json"))
    print(ccd_init)


def make_ccd_dirs(args: Namespace):
    os.chdir(args.calc_dir)
    gs, es = args.ccd_init.ground_structure, args.ccd_init.excited_structure
    e_to_g = es.interpolate(gs, nimages=args.e_to_g_div_ratios)
    g_to_e = gs.interpolate(es, nimages=args.g_to_e_div_ratios)

    e_dQs = [args.ccd_init.dQ * r for r in args.e_to_g_div_ratios]
    g_dQs = [args.ccd_init.dQ * (1.0 - r) for r in args.g_to_e_div_ratios]

    for state, ratios, structures, dQs in \
            [("excited", args.e_to_g_div_ratios, e_to_g, e_dQs),
             ("ground", args.g_to_e_div_ratios, g_to_e, g_dQs)]:

        if state == "ground":
            charge = args.ccd_init.ground_charge
            correction = args.ccd_init.ground_correction
        else:
            charge = args.ccd_init.excited_charge
            correction = args.ccd_init.excited_correction

        for ratio, structure, dQ in zip(ratios, structures, dQs):
            _make_dir(charge, state, ratio, structure, dQ, correction)


def _make_dir(charge, state, ratio, structure, dQ, correction):
    dir_ = Path(state) / f"disp_{ratio}"
    try:
        dir_.mkdir(parents=True, exist_ok=True)
        logger.info(f"Directory {dir_} was created.")

        structure.to(filename=dir_ / "POSCAR")
        (dir_ / "prior_info.yaml").write_text(
            yaml.dump({"charge": charge}), None)
        image_structure_info = ImageStructureInfo(dQ=dQ,
                                                  correction=correction,
                                                  correction_type="eFNV")
        image_structure_info.to_json_file(dir_ / "image_structure_info.json")

    except FileExistsError:
        logger.info(f"Directory {dir_} exists, so skip it.")


def make_ccd(args: Namespace):

    def _inner(dir_: Path):
        image_info: ImageStructureInfo = loadfn(dir_ / "image_structure_info.json")
        if image_info.energy is None:
            calc_results = get_calc_results(dir_, False)
            image_info.energy = calc_results.energy
        image_info.to_json_file(dir_ / "image_structure_info.json")
        return image_info

    ground_image_infos = parse_dirs(args.ground_dirs, _inner, verbose=True)
    excited_image_infos = parse_dirs(args.excited_dirs, _inner, verbose=True)

    ground_charge = args.ccd_init.ground_charge
    excited_charge = args.ccd_init.excited_charge

    if excited_charge - ground_charge == -1:
        ref = args.ccd_init.vbm
        majority_carrier_type = "p"
    else:
        ref = args.ccd_init.cbm
        majority_carrier_type = "n"

    for v in ground_image_infos:
        v.energy += ground_charge * ref
    for v in excited_image_infos:
        v.energy += excited_charge * ref

    ground_eg_image_infos = deepcopy(ground_image_infos)

    for v in ground_eg_image_infos:
        v.energy += args.ccd_init.cbm - args.ccd_init.vbm

    ccd = Ccd({"ground": ground_image_infos,
               f"excited + {majority_carrier_type}": excited_image_infos,
               "ground + p + n": ground_eg_image_infos},
              name=args.ccd_init.name)
    ccd.to_json_file()


def plot_ccd(args: Namespace):
    plotter = CcdPlotter(args.ccd, spline_deg=args.spline_deg)
    plotter.construct_plot()
    plotter.plt.savefig(args.fig_name)
    plotter.plt.show()


def plot_eigenvalues(args: Namespace):
    pass
#     qs, orb_infos = [], []
#     state = Path(args.dir).name
#
#     for d in glob(f'{args.dir}/disp_*'):
#         try:
#             orb_infos.append(loadfn(Path(d) / "band_edge_orbital_infos.json"))
#             disp_ratio = float(d.split("_")[-1])
#             qs.append(args.ccd.get_dQ_from_disp_ratio(state, disp_ratio))
#         except FileNotFoundError:
#             logger.info(f"band_edge_orbital_infos.json does not exist in {d}")
#             pass
#
#     pcr: PerfectBandEdgeState = args.perfect_band_edge_state
#     vbm, cbm = pcr.vbm_info.energy, pcr.cbm_info.energy
#     eigval_plotter = EigenvaluePlotter(orb_infos, qs, vbm, cbm)
#     eigval_plotter.construct_plot()
#     eigval_plotter.plt.savefig(args.fig_name)
#     eigval_plotter.plt.show()


