# -*- coding: utf-8 -*-
#  Copyright (c) 2022 Kumagai group.
import argparse
import sys
import warnings
from pathlib import Path

from monty.serialization import loadfn
from pymatgen.io.vasp.inputs import UnknownPotcarWarning

from dephon.version import __version__
from dephon.cli.main_function import make_ccd_init, make_ccd, \
    make_ccd_dirs, plot_ccd, plot_eigenvalues

warnings.simplefilter('ignore', UnknownPotcarWarning)


description = """Helper package to calculate carrier capture rates trapped by
point defects."""

epilog = f"Author: Yu Kumagai Version: {__version__}"


def parse_args_main(args):

    parser = argparse.ArgumentParser(description=description, epilog=epilog)
    subparsers = parser.add_subparsers()

    # -- make_ccd_init -----------------------------------
    parser_make_ccd_init = subparsers.add_parser(
        name="make_ccd_init",
        description="",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        aliases=['mci'])

    parser_make_ccd_init.add_argument(
        "-ed", "--excited_dir", type=Path, required=True,
        help="Directory for an excited state defect, e.g., Va_O1_0.")
    parser_make_ccd_init.add_argument(
        "-gd", "--ground_dir", type=Path, required=True,
        help="Directory for a ground state defect, e.g., Va_O1_1.")
    parser_make_ccd_init.set_defaults(func=make_ccd_init)

    # -- make_ccd_dirs -----------------------------------
    parser_add_ccd_dirs = subparsers.add_parser(
        name="make_ccd_dirs",
        description="",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        aliases=['mcd'])

    parser_add_ccd_dirs.add_argument(
        "--ccd_init", type=loadfn,
        default="ccd_init.json")
    parser_add_ccd_dirs.add_argument(
        "-egr", "--e_to_g_div_ratios", type=float, nargs="+",
        default=[-0.4, -0.2, 0.2, 0.4, 0.6, 0.8],
        help="Dividing ratios from excited state to ground state structures."
             "Thus, 1.0 means ground state structure")
    parser_add_ccd_dirs.add_argument(
        "-ger", "--g_to_e_div_ratios", type=float, nargs="+",
        default=[-0.4, -0.2, 0.2, 0.4, 0.6, 0.8],
        help="Dividing ratios from ground state to excited state structures.")
    parser_add_ccd_dirs.add_argument(
        "-d", "--calc_dir", type=Path, default=Path.cwd(),
        help="Directory where ground and excited directories are created.")
    parser_add_ccd_dirs.set_defaults(
        func=make_ccd_dirs)

    # -- make_ccd -----------------------------------
    parser_make_ccd = subparsers.add_parser(
        name="make_ccd",
        description="",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        aliases=['mc'])

    parser_make_ccd.add_argument(
        "--ccd_init", type=loadfn, default="ccd_init.json")
    parser_make_ccd.add_argument(
        "-g", "--ground_dirs", type=Path, nargs="+",
        help="Directories for ground directories.")
    parser_make_ccd.add_argument(
        "-e", "--excited_dirs", type=Path, nargs="+",
        help="Directories for excited directories.")
    parser_make_ccd.set_defaults(func=make_ccd)

    # -- plot_ccd -----------------------------------
    parser_plot_ccd = subparsers.add_parser(
        name="plot_ccd",
        description="",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        aliases=['pc'])

    parser_plot_ccd.add_argument(
        "--ccd", type=loadfn, default="ccd.json")
    parser_plot_ccd.add_argument(
        "--spline_deg", type=int, default=3)
    parser_plot_ccd.add_argument(
        "--fig_name", type=str, default="ccd.pdf")

    parser_plot_ccd.set_defaults(func=plot_ccd)

    # -- plot_eigenvalues -----------------------------------
    parser_plot_eigenvalues = subparsers.add_parser(
        name="plot_eigenvalues",
        description="",
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
    # ------------------------------------------------------------------------
    return parser.parse_args(args)


def main():
    args = parse_args_main(sys.argv[1:])
    args.func(args)


if __name__ == "__main__":
    main()




