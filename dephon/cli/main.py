# -*- coding: utf-8 -*-
#  Copyright (c) 2022 Kumagai group.
import argparse
import sys
import warnings
from pathlib import Path

from monty.serialization import loadfn
from pydefect.cli.main import add_sub_parser
from pymatgen.io.vasp.inputs import UnknownPotcarWarning

from dephon.cli.main_function import make_dephon_init, make_ccd, \
    make_ccd_dirs, plot_ccd, plot_eigenvalues, set_quadratic_fitting_q_range, \
    make_wswq_dirs, make_single_point_infos, make_single_ccd
from dephon.version import __version__

warnings.simplefilter('ignore', UnknownPotcarWarning)


description = """Helper package to calculate carrier capture rates trapped by
point defects."""

epilog = f"Author: Yu Kumagai Version: {__version__}"


def parse_args_main(args):

    parser = argparse.ArgumentParser(description=description, epilog=epilog)
    subparsers = parser.add_subparsers()

    unitcell_parser = add_sub_parser(argparse, name="unitcell")
    pbes_parser = add_sub_parser(argparse, name="perfect_band_edge_state")
    dirs = add_sub_parser(argparse, name="dirs")

    dephon_init = argparse.ArgumentParser(description="", add_help=False)
    dephon_init.add_argument(
            "--dephon_init", type=loadfn, default="dephon_init.json")

    # -- make_dephon_init -----------------------------------
    parser_make_dephon_init = subparsers.add_parser(
        name="make_dephon_init",
        description="Make dephon_init.json file from two directories with "
                    "pydefect files. When the excited state has one more "
                    "(less) charge state, n(p)-type is assumed.",
        parents=[unitcell_parser, pbes_parser],
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        aliases=['mdi'])

    parser_make_dephon_init.add_argument(
        "-fd", "--first_dir", type=Path, required=True,
        help="First directory considered for ccd, e.g., Va_O1_0.")
    parser_make_dephon_init.add_argument(
        "-sd", "--second_dir", type=Path, required=True,
        help="Second directory considered for ccd, e.g., Va_O1_1.")
    parser_make_dephon_init.set_defaults(func=make_dephon_init)

    # -- make_ccd_dirs -----------------------------------
    parser_add_ccd_dirs = subparsers.add_parser(
        name="make_ccd_dirs",
        description="Make directories to calculate configuration coordination "
                    "diagrams for ground and excited states.",
        parents=[dephon_init],
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        aliases=['mcd'])

    parser_add_ccd_dirs.add_argument(
        "-fsr", "--first_to_second_div_ratios", type=float, nargs="+",
        default=[-0.2, 0.0, 0.2, 0.4, 0.6, 0.8],
        help="Dividing ratios from first to second charge state structures.")
    parser_add_ccd_dirs.add_argument(
        "-sfr", "--second_to_first_div_ratios", type=float, nargs="+",
        default=[-0.2, 0.0, 0.2, 0.4, 0.6, 0.8],
        help="Dividing ratios from second to first charge state structures.")
    parser_add_ccd_dirs.add_argument(
        "-d", "--calc_dir", type=Path, default=Path.cwd(),
        help="Directory where directories are created.")
    parser_add_ccd_dirs.set_defaults(
        func=make_ccd_dirs)

    # -- make_single_point_infos -----------------------------------
    parser_make_single_point_infos = subparsers.add_parser(
        name="make_single_point_infos",
        description="Make single_point_info.json at each directory. "
                    "Before running this command, calc_results.json and "
                    "band_edge_states.json files need to be created using "
                    "pydefect.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        parents=[dirs],
        aliases=['mspi'])

    parser_make_single_point_infos.set_defaults(func=make_single_point_infos)

    # -- make_single_ccd -----------------------------------
    parser_make_single_ccd = subparsers.add_parser(
        name="make_single_ccd",
        description="Make single_point_info.json at each directory. "
                    "Before running this command, calc_results.json and "
                    "band_edge_states.json files need to be created using "
                    "pydefect.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        parents=[dirs],
        aliases=['msc'])

    parser_make_single_ccd.set_defaults(func=make_single_ccd)

    # -- make_ccd -----------------------------------
    parser_make_ccd = subparsers.add_parser(
        name="make_ccd",
        description="Make ccd.json file",
        parents=[dephon_init],
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        aliases=['mc'])

    parser_make_ccd.add_argument(
        "-g", "--ground_ccd", type=loadfn,
        help="single_ccd.json file corresponding to a ground state.")
    parser_make_ccd.add_argument(
        "-e", "--excited_ccd", type=loadfn,
        help="single_ccd.json file corresponding to an excited state.")
    parser_make_ccd.set_defaults(func=make_ccd)

    # -- set_fitting_q_range -----------------------------------
    parser_set_fitting_q_range = subparsers.add_parser(
        name="set_fitting_q_range",
        description="Set the fitting range for the quadratic potential surface",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        aliases=['sfr'])

    parser_set_fitting_q_range.add_argument(
        "--ccd", type=loadfn, default="ccd.json")
    parser_set_fitting_q_range.add_argument(
        "--single_ccd_name", type=str, required=True)
    parser_set_fitting_q_range.add_argument(
        "--q_range", type=float, nargs="+")
    parser_set_fitting_q_range.set_defaults(func=set_quadratic_fitting_q_range)

    # -- plot_ccd -----------------------------------
    parser_plot_ccd = subparsers.add_parser(
        name="plot_ccd",
        description="Plot cc diagram from ccd.json file",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        aliases=['pc'])

    parser_plot_ccd.add_argument(
        "--ccd", type=loadfn, default="ccd.json")
    parser_plot_ccd.add_argument(
        "--fig_name", type=str, default="ccd.pdf")
    parser_plot_ccd.add_argument(
        "--q_range", type=float, nargs="+")

    parser_plot_ccd.set_defaults(func=plot_ccd)

    # -- plot_eigenvalues -----------------------------------
    parser_plot_eigenvalues = subparsers.add_parser(
        name="plot_eigenvalues",
        parents=[dephon_init, dirs],
        description="Plot eigenvalues as function of displacement ratio. "
                    "band_edge_orbital_infos.json is needed at each directory.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        aliases=['pe'])

    parser_plot_eigenvalues.set_defaults(func=plot_eigenvalues)

    # -- make_wswq_dirs -----------------------------------
    parser_make_wswq_dirs = subparsers.add_parser(
        name="make_wswq_dirs",
        description="Make directories for calculating WSWQ files.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        aliases=['mwd'])

    parser_make_wswq_dirs.add_argument(
        "--ccd_init", type=loadfn, default="ccd_init.json")
    parser_make_wswq_dirs.add_argument(
        "--ground_dirs", type=Path, nargs="+", default=[])
    parser_make_wswq_dirs.add_argument(
        "--excited_dirs", type=Path, nargs="+", default=[])

    parser_make_wswq_dirs.set_defaults(func=make_wswq_dirs)
    # ------------------------------------------------------------------------
    return parser.parse_args(args)


def main():
    args = parse_args_main(sys.argv[1:])
    args.func(args)


if __name__ == "__main__":
    main()




