# -*- coding: utf-8 -*-
#  Copyright (c) 2022 Kumagai group.
import os
from argparse import Namespace
from pathlib import Path

import yaml
from monty.serialization import loadfn
from pydefect.analyzer.band_edge_states import BandEdgeStates, \
    BandEdgeOrbitalInfos
from pydefect.analyzer.calc_results import CalcResults
from pydefect.analyzer.defect_energy import DefectEnergyInfo
from pydefect.analyzer.defect_structure_info import DefectStructureInfo
from pydefect.cli.main_functions import get_calc_results
from pydefect.cli.main_tools import parse_dirs
from pymatgen.electronic_structure.core import Spin
from vise.input_set.incar import ViseIncar
from vise.input_set.prior_info import PriorInfo
from vise.util.file_transfer import FileLink
from vise.util.logger import get_logger

from dephon.config_coord import SinglePointInfo, CcdPlotter, \
    SingleCcd, SingleCcdId
from dephon.corrections import DephonCorrection
from dephon.dephon_init import DephonInit, MinimumPointInfo, NearEdgeState
from dephon.enum import CorrectionType
from dephon.make_config_coord import MakeCcd
from dephon.make_ele_phon_coupling import MakeInitialEPCoupling
from dephon.plot_eigenvalues import DephonEigenvaluePlotter
from dephon.util import spin_to_idx

logger = get_logger(__name__)


def make_near_edge_states(band_edge_orbital_infos: BandEdgeOrbitalInfos,
                          spin: Spin,
                          edge_energy: float,
                          threshold: float = 0.1):
    result = []
    info_by_spin = band_edge_orbital_infos.orbital_infos[spin_to_idx(spin)]
    for k_idx, (orb_info_by_kpt, k_coords, k_weight) in \
        enumerate(zip(info_by_spin,
                      band_edge_orbital_infos.kpt_coords,
                      band_edge_orbital_infos.kpt_weights)):
        k_idx_from_1 = k_idx + 1
        for rel_idx, info_by_band in enumerate(orb_info_by_kpt):
            e_from_band_edge = abs(info_by_band.energy - edge_energy)
            if e_from_band_edge > threshold:
                continue
            band_idx = rel_idx + band_edge_orbital_infos.lowest_band_index
            result.append(NearEdgeState(band_idx,
                                        k_coords,
                                        k_weight,
                                        k_idx_from_1,
                                        info_by_band.energy,
                                        info_by_band.occupation))
    return result


def make_min_point_info_from_dir(_dir: Path):
    energy_info = DefectEnergyInfo.from_yaml(_dir / "defect_energy_info.yaml")
    calc_results: CalcResults = loadfn(_dir / "calc_results.json")
    defect_structure_info: DefectStructureInfo \
        = loadfn(_dir / "defect_structure_info.json")
    band_edge_states: BandEdgeStates = loadfn(_dir / "band_edge_states.json")

    band_edge_orbital_infos: BandEdgeOrbitalInfos = loadfn(_dir / "band_edge_orbital_infos.json")
    localized_orbitals, valence_bands, conduction_bands = [], [], []
    for state, spin in zip(band_edge_states.states, [Spin.up, Spin.down]):
        localized_orbitals.append(state.localized_orbitals)
        valence_bands.append(make_near_edge_states(band_edge_orbital_infos,
                                                   spin,
                                                   state.vbm_info.energy))
        conduction_bands.append(make_near_edge_states(band_edge_orbital_infos,
                                                      spin,
                                                      state.cbm_info.energy))
    min_point_info = MinimumPointInfo(
        charge=energy_info.charge,
        structure=calc_results.structure,
        energy=energy_info.defect_energy.formation_energy,
        correction_energy=energy_info.defect_energy.total_correction,
        magnetization=calc_results.magnetization,
        localized_orbitals=localized_orbitals,
        initial_site_symmetry=defect_structure_info.initial_site_sym,
        final_site_symmetry=defect_structure_info.final_site_sym,
        parsed_dir=str(_dir.absolute()),
        valence_bands=valence_bands,
        conduction_bands=conduction_bands)

    return min_point_info, energy_info.name


def make_dephon_init(args: Namespace):
    state1, name1 = make_min_point_info_from_dir(args.first_dir)
    state2, name2 = make_min_point_info_from_dir(args.second_dir)

    if name1 != name2:
        logger.warning("The names of ground and excited states are "
                       f"{name1} and {name2}. Here, {name1} is used.")

    concentration = args.effective_mass.concentrations[0]
    ave_hole_mass = args.effective_mass.average_mass("p", concentration)
    ave_electron_mass = args.effective_mass.average_mass("n", concentration)
    dephon_init = DephonInit(defect_name=name1,
                             states=[state1, state2],
                             vbm=args.unitcell.vbm,
                             cbm=args.unitcell.cbm,
                             supercell_vbm=args.p_state.vbm_info.energy,
                             supercell_cbm=args.p_state.cbm_info.energy,
                             ave_hole_mass=ave_hole_mass,
                             ave_electron_mass=ave_electron_mass,
                             ave_static_diele_const=args.unitcell.ave_ele_diele)

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
            id_ = SingleCcdId(name)
            SingleCcd(id_, charge=initial_charge).to_json_file(single_ccd)


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


def update_single_point_infos(args: Namespace):
    def _inner(dir_: Path):
        calc_results = get_calc_results(dir_, False)
        band_edge_states: BandEdgeStates = loadfn(dir_ / "band_edge_states.json")
        correction = DephonCorrection.from_yaml(dir_ / "dephon_correction.yaml")

        sp_info: SinglePointInfo = loadfn(dir_ / "single_point_info.json")
        sp_info.corrected_energy = calc_results.energy + correction.energy
        sp_info.magnetization = calc_results.magnetization
        sp_info.localized_orbitals = [s.localized_orbitals for s in band_edge_states.states]
        sp_info.is_shallow = band_edge_states.is_shallow
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
    ccd = MakeCcd(args.ground_ccd, args.excited_ccd, args.dephon_init).ccd
    print(ccd)
    ccd.to_json_file()


def set_quadratic_fitting_q_range(args: Namespace):
    single_ccd: SingleCcd = args.ccd.single_ccd(args.single_ccd_name)
    single_ccd.set_quadratic_fitting_range(args.q_range)
    print(args.ccd)
    args.ccd.to_json_file()


def plot_ccd(args: Namespace):
    plotter = CcdPlotter(args.ccd, q_range=args.q_range)
    plotter.construct_plot()
    plotter.plt.savefig(args.fig_name)
    plotter.plt.show()


def plot_eigenvalues(args: Namespace):
    disp_ratios, orb_infos = [], []

    for d in args.dirs:
        try:
            orb_info = loadfn(Path(d) / "band_edge_orbital_infos.json")
        except FileNotFoundError:
            logger.info(f"band_edge_orbital_infos.json does not exist in {d}")
            continue
        try:
            single_point_info = loadfn(Path(d) / "single_point_info.json")
        except FileNotFoundError:
            logger.info(f"single_point_info.json does not exist in {d}")
            continue
        orb_infos.append(orb_info)
        disp_ratios.append(single_point_info.disp_ratio)

    vbm, cbm = args.dephon_init.supercell_vbm, args.dephon_init.supercell_cbm
    eigval_plotter = DephonEigenvaluePlotter(orb_infos, disp_ratios, vbm, cbm)
    eigval_plotter.construct_plot()
    eigval_plotter.plt.savefig("dephon_eigenvalues.pdf")
    eigval_plotter.plt.show()


def make_wswq_dirs(args: Namespace):
    for dir_ in args.dirs:
        _make_wswq_dir(dir_, args.dephon_init)


def _make_wswq_dir(dir_, dephon_init: DephonInit):
    wswq_dir = (dir_ / "wswq")
    if wswq_dir.exists():
        logger.info(f"Directory {wswq_dir} exists, so skip creating it.")
        return

    charge = PriorInfo.load_yaml(dir_ / "prior_info.yaml").charge
    original_dir = dephon_init.min_info_from_charge(charge).dir_path

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


def make_initial_e_p_coupling(args: Namespace):
    make_init = MakeInitialEPCoupling(args.dephon_init,
                                      args.ccd,
                                      args.captured_carrier,
                                      args.charge_for_e_p_coupling)
    e_p_coupling = make_init.make()
    print(e_p_coupling)
    e_p_coupling.to_json_file()
#
#
# def update_e_p_coupling(args: Namespace):
#     add_inner_products(args.)