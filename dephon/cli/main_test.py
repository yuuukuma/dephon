# -*- coding: utf-8 -*-
#  Copyright (c) 2022 Kumagai group.
from argparse import Namespace
from pathlib import Path

from dephon.cli.main import parse_args_main


def test_main_make_ccd_init_and_dirs():
    parsed_args = parse_args_main(["ci",
                                   "-id", "Va_O1_0",
                                   "-fd", "Va_O1_1"])
    expected = Namespace(
        initial_dir=Path("Va_O1_0"),
        final_dir=Path("Va_O1_1"),
        i_to_f_div_ratios=[-0.4, -0.2, 0.2, 0.4, 0.6, 0.8],
        f_to_i_div_ratios=[-0.4, -0.2, 0.2, 0.4, 0.6, 0.8],
        func=parsed_args.func)
    assert parsed_args == expected


def test_main_make_ccd_init_and_dirs_w_args():
    parsed_args = parse_args_main(["ci",
                                   "-id", "Va_O1_3",
                                   "-fd", "Va_O1_4",
                                   "-ifr", "0.1",
                                   "-fir", "-0.1"])
    expected = Namespace(
        initial_dir=Path("Va_O1_3"),
        final_dir=Path("Va_O1_4"),
        i_to_f_div_ratios=[0.1],
        f_to_i_div_ratios=[-0.1],
        func=parsed_args.func)
    assert parsed_args == expected

