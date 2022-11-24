# -*- coding: utf-8 -*-
#  Copyright (c) 2022 Kumagai group.
import argparse
import sys
import warnings
from pathlib import Path

from monty.serialization import loadfn
from pymatgen.io.vasp.inputs import UnknownPotcarWarning

from dephon.version import __version__
from dephon.cli.main_function import make_ccd_init_and_dirs, make_ccd, \
    add_ccd_dirs

warnings.simplefilter('ignore', UnknownPotcarWarning)


description = """Helper package to calculate carrier capture rates trapped by
point defects."""

epilog = f"Author: Yu Kumagai Version: {__version__}"


def parse_args_main(args):

    parser = argparse.ArgumentParser(description=description, epilog=epilog)
    subparsers = parser.add_subparsers()

    # -- make_ccd_init_and_dirs -----------------------------------
    parser_make_ccd_init = subparsers.add_parser(
        name="make_ccd_init_and_dirs",
        description="",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        aliases=['ci'])

    parser_make_ccd_init.add_argument(
        "-ed", "--excited_dir", type=Path, required=True,
        help="Directory for an excited state defect.")
    parser_make_ccd_init.add_argument(
        "-gd", "--ground_dir", type=Path, required=True,
        help="Directory for a ground state defect.")
    parser_make_ccd_init.add_argument(
        "-egr", "--e_to_g_div_ratios", type=float, nargs="+",
        default=[-0.4, -0.2, 0.2, 0.4, 0.6, 0.8],
        help="Dividing ratios from excited state to ground state structures.")
    parser_make_ccd_init.add_argument(
        "-ger", "--g_to_e_div_ratios", type=float, nargs="+",
        default=[-0.4, -0.2, 0.2, 0.4, 0.6, 0.8],
        help="Dividing ratios from ground state to excited state structures.")
    parser_make_ccd_init.set_defaults(
        func=make_ccd_init_and_dirs)

    # -- add_ccd_dirs -----------------------------------
    parser_add_ccd_dirs = subparsers.add_parser(
        name="add_ccd_dirs",
        description="",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        aliases=['ac'])

    parser_add_ccd_dirs.add_argument(
        "--ccd_init", type=loadfn,
        default="ccd_init.json")
    parser_add_ccd_dirs.add_argument(
        "-egr", "--e_to_g_div_ratios", type=float, nargs="+",
        default=[],
        help="Dividing ratios from excited state to ground state structures.")
    parser_add_ccd_dirs.add_argument(
        "-ger", "--g_to_e_div_ratios", type=float, nargs="+",
        default=[],
        help="Dividing ratios from ground state to excited state structures.")
    parser_add_ccd_dirs.add_argument(
        "-d", "--calc_dir", type=Path, nargs="+",
        default=Path.cwd(),
        help="Directory where ground and excited directories are created.")
    parser_add_ccd_dirs.set_defaults(
        func=add_ccd_dirs)

    # -- make_ccd -----------------------------------
    parser_make_ccd = subparsers.add_parser(
        name="make_ccd",
        description="",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        aliases=['c'])

    parser_make_ccd.add_argument(
        "--ccd_init", type=loadfn, default="ccd_init.json")
    parser_make_ccd.set_defaults(func=make_ccd)

    # ------------------------------------------------------------------------
    return parser.parse_args(args)


def main():
    args = parse_args_main(sys.argv[1:])
    args.func(args)


if __name__ == "__main__":
    main()




