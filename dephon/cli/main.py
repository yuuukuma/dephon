# -*- coding: utf-8 -*-
#  Copyright (c) 2022 Kumagai group.
import argparse
import sys
import warnings
from pathlib import Path

from pymatgen.io.vasp.inputs import UnknownPotcarWarning

from dephon import __version__
from dephon.cli.main_function import make_ccd_init_and_dirs

warnings.simplefilter('ignore', UnknownPotcarWarning)


description = """Helper package to calculate carrier capture rates trapped by
point defects."""

epilog = f"Author: Yu Kumagai Version: {__version__}"


def parse_args_main(args):

    parser = argparse.ArgumentParser(description=description, epilog=epilog)
    subparsers = parser.add_subparsers()

    # -- make_ccd_init_and_dirs -----------------------------------
    parser_make_cpp_init = subparsers.add_parser(
        name="make_ccd_init_and_dirs",
        description="",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        aliases=['ci'])

    parser_make_cpp_init.add_argument(
        "-id", "--initial_dir", type=Path, required=True,
        help="Directory for an initial excited state defect.")
    parser_make_cpp_init.add_argument(
        "-fd", "--final_dir", type=Path, required=True,
        help="Directory for a final ground state defect.")
    parser_make_cpp_init.add_argument(
        "-ifr", "--i_to_f_div_ratios", type=float, nargs="+",
        default=[-0.4, -0.2, 0.2, 0.4, 0.6, 0.8],
        help="Dividing ratios from initial (excited state) to final "
             "(ground state) structures.")
    parser_make_cpp_init.add_argument(
        "-fir", "--f_to_i_div_ratios", type=float, nargs="+",
        default=[-0.4, -0.2, 0.2, 0.4, 0.6, 0.8],
        help="Dividing ratios from final (ground state) to initial "
             "(excited state) structures.")
    parser_make_cpp_init.set_defaults(
        func=make_ccd_init_and_dirs)
    # ------------------------------------------------------------------------
    return parser.parse_args(args)


def main():
    args = parse_args_main(sys.argv[1:])
    args.func(args)


if __name__ == "__main__":
    main()




