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
    make_ccd_dirs, plot_ccd, plot_eigenvalues, set_fitting_q_range, \
    make_wswq_dirs
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
        default=[-0.4, -0.2, 0.0, 0.2, 0.4, 0.6, 0.8],
        help="Dividing ratios from first to second charge state structures.")
    parser_add_ccd_dirs.add_argument(
        "-sfr", "--second_to_first_div_ratios", type=float, nargs="+",
        default=[-0.4, -0.2, 0.0, 0.2, 0.4, 0.6, 0.8],
        help="Dividing ratios from second to first charge state structures.")
    parser_add_ccd_dirs.add_argument(
        "-d", "--calc_dir", type=Path, default=Path.cwd(),
        help="Directory where ground and excited directories are created.")
    parser_add_ccd_dirs.set_defaults(
        func=make_ccd_dirs)

    # -- make_ccd -----------------------------------
    parser_make_ccd = subparsers.add_parser(
        name="make_ccd",
        description="Make ccd.json file from calculated directories."
                    "Before running this command, one needs to create "
                    "calc_results.json and band_edge_states.json files"
                    "at each directory.",
        parents=[dephon_init],
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        aliases=['mc'])

    parser_make_ccd.add_argument(
        "-g", "--ground_dirs", type=Path, nargs="+",
        help="Directories for ground directories.")
    parser_make_ccd.add_argument(
        "-e", "--excited_dirs", type=Path, nargs="+",
        help="Directories for excited directories.")
    parser_make_ccd.add_argument(
        "-s", "--skip_shallow", action="store_true",
        help="Set when skip shallow states.")
    parser_make_ccd.set_defaults(func=make_ccd)

    # -- set_fitting_q_range -----------------------------------
    parser_set_fitting_q_range = subparsers.add_parser(
        name="set_fitting_q_range",
        description="Set the fitting range for the quadratic potential surface",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        aliases=['sfq'])

    parser_set_fitting_q_range.add_argument(
        "--ccd", type=loadfn, default="ccd.json")
    parser_set_fitting_q_range.add_argument(
        "--image_name", type=str, required=True)
    parser_set_fitting_q_range.add_argument(
        "--q_min", type=float)
    parser_set_fitting_q_range.add_argument(
        "--q_max", type=float)
    parser_set_fitting_q_range.set_defaults(func=set_fitting_q_range)

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

    parser_plot_ccd.set_defaults(func=plot_ccd)

    # -- plot_eigenvalues -----------------------------------
    parser_plot_eigenvalues = subparsers.add_parser(
        name="plot_eigenvalues",
        parents=[dephon_init],
        description="Plot eigenvalues as function of Q for each state.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        aliases=['pe'])

    parser_plot_eigenvalues.add_argument(
        "--ccd", type=loadfn, default="ccd.json")
    parser_plot_eigenvalues.add_argument(
        "--dir", type=Path,
        help="Set ground or excited directory path. E.g. XXX/YY/ground")
    parser_plot_eigenvalues.add_argument(
        "--fig_name", type=str, default="eigenvalues.pdf")
    parser_plot_eigenvalues.add_argument(
        "-pbes", "--perfect_band_edge_state", type=loadfn,
        help="Path to the perfect_band_edge_state.json.")
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




