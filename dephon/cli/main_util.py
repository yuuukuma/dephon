# -*- coding: utf-8 -*-
#  Copyright (c) 2022 Kumagai group.
import argparse
import sys
import warnings
from pathlib import Path

from monty.serialization import loadfn
from pymatgen.io.vasp.inputs import UnknownPotcarWarning

from dephon.cli.main_util_function import reduce_wswq_auto
from dephon.version import __version__

warnings.simplefilter('ignore', UnknownPotcarWarning)


description = """Helper package for calculating the non-radiative carrier 
capture rates trapped by point defects."""

epilog = f"Author: Yu Kumagai Version: {__version__}"


def parse_args_main(args):

    parser = argparse.ArgumentParser(description=description, epilog=epilog)
    subparsers = parser.add_subparsers()

    dephon_init = argparse.ArgumentParser(description="", add_help=False)
    dephon_init.add_argument(
            "--dephon_init", type=loadfn, default="dephon_init.json")

    ccd = argparse.ArgumentParser(description="", add_help=False)
    ccd.add_argument(
        "--ccd", type=loadfn, default="ccd.json")

    # -- reduce_wswq_auto -----------------------------------
    parser_reduce_wswq_auto = subparsers.add_parser(
        name="reduce_wswq_auto",
        description="",
        parents=[dephon_init],
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        aliases=['rwa'])

    parser_reduce_wswq_auto.add_argument(
        "-s", "--single_ccd", type=loadfn, required=True,
        help="single_ccd.json file.")
    parser_reduce_wswq_auto.add_argument(
        "-w", "--wswq", type=Path, required=True,
        help="")
    parser_reduce_wswq_auto.add_argument(
        "-b", "--band_indices", type=float, nargs="+",
        help="")

    parser_reduce_wswq_auto.set_defaults(func=reduce_wswq_auto)
    # ------------------------------------------------------------------------
    return parser.parse_args(args)


def main():
    args = parse_args_main(sys.argv[1:])
    args.func(args)


if __name__ == "__main__":
    main()




