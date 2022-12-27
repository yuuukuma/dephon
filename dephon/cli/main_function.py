# -*- coding: utf-8 -*-
#  Copyright (c) 2022 Kumagai group.
import os
from argparse import Namespace
from copy import deepcopy
from glob import glob
from pathlib import Path

import yaml
from monty.serialization import loadfn
from pydefect.analyzer.band_edge_states import BandEdgeState
from pydefect.analyzer.calc_results import CalcResults
from pydefect.analyzer.defect_energy import DefectEnergyInfo
from pydefect.analyzer.defect_structure_info import DefectStructureInfo
from pydefect.cli.main_functions import get_calc_results
from pydefect.cli.main_tools import parse_dirs
from vise.input_set.incar import ViseIncar
from vise.util.file_transfer import FileLink
from vise.util.logger import get_logger

from dephon.config_coord import SinglePointInfo, Ccd, CcdPlotter, \
    SingleCcd
from dephon.corrections import DephonCorrection
from dephon.dephon_init import DephonInit, MinimumPointInfo
from dephon.enum import CorrectionType
from dephon.plot_eigenvalues import EigenvaluePlotter

logger = get_logger(__name__)


def _make_min_point_info_from_dir(_dir: Path):
    energy_info = DefectEnergyInfo.from_yaml(_dir / "defect_energy_info.yaml")
    calc_results: CalcResults = loadfn(_dir / "calc_results.json")
    defect_structure_info: DefectStructureInfo \
        = loadfn(_dir / "defect_structure_info.json")
    min_point_info = MinimumPointInfo(
        charge=energy_info.charge,
        structure=calc_results.structure,
        energy=energy_info.defect_energy.formation_energy,
        correction_energy=energy_info.defect_energy.total_correction,
        initial_site_symmetry=defect_structure_info.initial_site_sym,
        final_site_symmetry=defect_structure_info.final_site_sym,
        parsed_dir=str(_dir.absolute()))
    return min_point_info, energy_info.name


def make_dephon_init(args: Namespace):
    state1, name1 = _make_min_point_info_from_dir(args.first_dir)
    state2, name2 = _make_min_point_info_from_dir(args.second_dir)

    if name1 != name2:
        logger.warning("The names of ground and excited states are "
                       f"{name1} and {name2}. Here, {name1} is used.")
    dephon_init = DephonInit(defect_name=name1,
                             states=[state1, state2],
                             vbm=args.unitcell.vbm,
                             cbm=args.unitcell.cbm,
                             supercell_vbm=args.p_state.vbm_info.energy,
                             supercell_cbm=args.p_state.cbm_info.energy)

    path = Path(f"cc/{dephon_init.defect_name}_{state1.charge}_{state2.charge}")
    if path.exists() is False:
        path.mkdir(parents=True)

    json_file = path / "dephon_init.json"
    if json_file.exists():
        logger.info(f"{json_file} exists. Remove it first to recreate it.")
        return

    dephon_init.to_json_file(json_file)
    print(dephon_init)


def make_ccd_dirs(args: Namespace):
    os.chdir(args.calc_dir)
    d_init = args.dephon_init
    s1 = d_init.states[0].structure
    s2 = d_init.states[1].structure
    s1_to_s2 = s1.interpolate(s2, nimages=args.first_to_second_div_ratios)
    s2_to_s1 = s2.interpolate(s1, nimages=args.second_to_first_div_ratios)
    initial_charges = [d_init.states[0].charge, d_init.states[1].charge]
    final_charges = [d_init.states[1].charge, d_init.states[0].charge]
    correction_energies = [d_init.states[0].correction_energy,
                           d_init.states[1].correction_energy]

    for ratios, structures, initial_charge, final_charge, corr in \
            zip([args.first_to_second_div_ratios,
                 args.second_to_first_div_ratios],
                [s1_to_s2, s2_to_s1],
                initial_charges, final_charges, correction_energies):
        dQs = [args.dephon_init.dQ * r for r in ratios]
        name = f"from_{initial_charge}_to_{final_charge}"

        for ratio, structure, dQ in zip(ratios, structures, dQs):
            _make_ccd_dir(initial_charge, name, ratio, structure, dQ, corr)

        single_ccd = Path(name) / "single_ccd.json"
        if single_ccd.exists() is False:
            SingleCcd(name=name, charge=initial_charge).to_json_file(single_ccd)


def _make_ccd_dir(charge, dirname, ratio, structure, dQ, correction):
    dir_ = Path(dirname) / f"disp_{ratio}"
    try:
        dir_.mkdir(parents=True)
        logger.info(f"Directory {dir_} was created.")

        structure.to(filename=str(dir_ / "POSCAR"))
        (dir_ / "prior_info.yaml").write_text(
            yaml.dump({"charge": charge}), None)
        single_point_info = SinglePointInfo(dQ=dQ, disp_ratio=ratio)
        single_point_info.to_json_file(dir_ / "single_point_info.json")

        correction = DephonCorrection(correction, CorrectionType.extended_FNV)
        correction.to_yaml_file(dir_ / "dephon_correction.yaml")

    except FileExistsError:
        logger.info(f"Directory {dir_} exists, so skip it.")


def make_single_point_infos(args: Namespace):
    def _inner(dir_: Path):
        calc_results = get_calc_results(dir_, False)
        band_edge_state: BandEdgeState = loadfn(dir_ / "band_edge_states.json")
        correction = DephonCorrection.from_yaml(dir_ / "dephon_correction.yaml")

        sp_info: SinglePointInfo = loadfn(dir_ / "single_point_info.json")
        sp_info.corrected_energy = calc_results.energy + correction.energy
        sp_info.is_shallow = band_edge_state.is_shallow
        sp_info.correction_method = correction.correction_type
        sp_info.to_json_file(dir_ / "single_point_info.json")

    parse_dirs(args.dirs, _inner, verbose=True)


def make_single_ccd(args: Namespace):
    def _inner(dir_: Path):
        return loadfn(dir_ / "single_point_info.json")

    single_ccd: SingleCcd = loadfn("single_ccd.json")
    single_ccd.point_infos = parse_dirs(args.dirs, _inner, verbose=True)
    single_ccd.to_json_file("single_ccd.json")


def make_ccd(args: Namespace):

    def _inner(dir_: Path):
        image_info: SinglePointInfo = loadfn(dir_ / "image_structure_info.json")
        if image_info.bare_energy is None:
            calc_results = get_calc_results(dir_, False)
            image_info.bare_energy = calc_results.energy
        try:
            band_edge_state: BandEdgeState = loadfn(dir_ / "band_edge_states.json")
            image_info.is_shallow = band_edge_state.is_shallow
            if band_edge_state.is_shallow:
                logger.info(f"{dir_} has shallow carriers.")
        except FileNotFoundError:
            logger.warning("To judge if the states are shallow or not,"
                           "we need band_edge_states.json.")
        image_info.to_json_file(dir_ / "image_structure_info.json")
        return image_info

    ground_image_infos = parse_dirs(args.ground_dirs, _inner, verbose=True)
    excited_image_infos = parse_dirs(args.excited_dirs, _inner, verbose=True)
    ground_eg_image_infos = deepcopy(ground_image_infos)

    _add_carrier_energies(args.dephon_init, ground_image_infos,
                          excited_image_infos, ground_eg_image_infos)

    excited_name = f"excited + {args.dephon_init.semiconductor_type}"
    ground = SingleCcd("ground", ground_image_infos)
    excited = SingleCcd(excited_name, excited_image_infos)
    ground_eg = SingleCcd("ground + p + n", ground_eg_image_infos)

    ccd = Ccd(defect_name="test",
              correction_energy_type=CorrectionType.extended_FNV,
              image_infos_list=[ground, excited, ground_eg])
    ccd.to_json_file()


def _add_carrier_energies(dephon_init, ground_image_infos, excited_image_infos,
                          ground_eg_image_infos):
    fermi_level_from_vbm = dephon_init.delta_EF
    for v in ground_image_infos:
        v.bare_energy += dephon_init.single_ccd.charge * fermi_level_from_vbm
    for v in excited_image_infos:
        v.bare_energy += dephon_init.excited_state.charge * fermi_level_from_vbm
    for v in ground_eg_image_infos:
        v.bare_energy += dephon_init.band_gap


def set_fitting_q_range(args: Namespace):
    ccd: Ccd = args.ccd
    imag = ccd.single_ccd(args.image_name)
    imag.set_q_range(args.q_min, args.q_max)
    ccd.to_json_file()
    print(imag)


def plot_ccd(args: Namespace):
    plotter = CcdPlotter(args.ccd)
    plotter.construct_plot()
    plotter.plt.savefig(args.fig_name)
    plotter.plt.show()


def plot_eigenvalues(args: Namespace):
    qs, orb_infos = [], []

    for d in glob(f'{args.dir}/disp_*'):
        try:
            orb_infos.append(loadfn(Path(d) / "band_edge_orbital_infos.json"))
            imag_info: SinglePointInfo = loadfn(Path(d) / "image_structure_info.json")
            qs.append(imag_info.dQ)
        except FileNotFoundError:
            logger.info(f"band_edge_orbital_infos.json does not exist in {d}")
            pass

    vbm, cbm = args.dephon_init.vbm, args.dephon_init.cbm
    eigval_plotter = EigenvaluePlotter(orb_infos, qs, vbm, cbm)
    eigval_plotter.construct_plot()
    eigval_plotter.plt.savefig(args.fig_name)
    eigval_plotter.plt.show()


def make_wswq_dirs(args: Namespace):
    for dir_ in args.ground_dirs:
        _make_wswq_dir(dir_, args.dephon_init.single_ccd.dir_path)
    for dir_ in args.excited_dirs:
        _make_wswq_dir(dir_, args.dephon_init.excited_state.dir_path)


def _make_wswq_dir(dir_, original_dir):
    wswq_dir = (dir_ / "wswq")
    if wswq_dir.exists():
        logger.info(f"Directory {wswq_dir} exists, so skip creating it.")
        return

    wswq_dir.mkdir()
    logger.info(f"Directory {wswq_dir} was created.")

    for f_name in ["KPOINTS", "POSCAR", "POTCAR"]:
        FileLink((dir_/f_name).absolute()).transfer(wswq_dir)

    incar = ViseIncar.from_file(dir_/"INCAR")
    incar.update({"ALGO": "None", "LWSWQ": True, "NELM": 1, "LWAVE": False})
    incar.pop("LORBIT", None)
    incar.write_file(Path(wswq_dir/"INCAR"))

    os.symlink((dir_/"WAVECAR").absolute(), (wswq_dir/"WAVECAR.qqq"))

    os.symlink((original_dir/"WAVECAR").absolute(), (wswq_dir/"WAVECAR"))
