# -*- coding: utf-8 -*-
#  Copyright (c) 2022 Kumagai group.
from argparse import Namespace
from pathlib import Path

from dephon.cli.main import parse_args_main
from dephon.config_coord import CcdInit, Ccd


def test_main_make_ccd_init_and_dirs():
    parsed_args = parse_args_main(["ci",
                                   "-ed", "Va_O1_0",
                                   "-gd", "Va_O1_1"])
    expected = Namespace(
        excited_dir=Path("Va_O1_0"),
        ground_dir=Path("Va_O1_1"),
        e_to_g_div_ratios=[-0.4, -0.2, 0.2, 0.4, 0.6, 0.8],
        g_to_e_div_ratios=[-0.4, -0.2, 0.2, 0.4, 0.6, 0.8],
        func=parsed_args.func)
    assert parsed_args == expected


def test_main_make_ccd_init_and_dirs_w_args():
    parsed_args = parse_args_main(["ci",
                                   "-ed", "Va_O1_3",
                                   "-gd", "Va_O1_4",
                                   "-egr", "0.1",
                                   "-ger", "-0.1"])
    expected = Namespace(
        excited_dir=Path("Va_O1_3"),
        ground_dir=Path("Va_O1_4"),
        e_to_g_div_ratios=[0.1],
        g_to_e_div_ratios=[-0.1],
        func=parsed_args.func)
    assert parsed_args == expected


def test_main_make_ccd_wo_args(mocker):
    mock_ccd_init = mocker.Mock(spec=CcdInit, autospec=True)

    def side_effect(filename):
        if filename == "ccd_init.json":
            return mock_ccd_init
        else:
            raise ValueError

    mocker.patch("dephon.cli.main.loadfn", side_effect=side_effect)

    parsed_args = parse_args_main(["c"])
    expected = Namespace(
        ccd_init=mock_ccd_init,
        func=parsed_args.func)
    assert parsed_args == expected


def test_main_plot_ccd_wo_args(mocker):
    mock_ccd = mocker.Mock(spec=Ccd, autospec=True)

    def side_effect(filename):
        if filename == "ccd.json":
            return mock_ccd
        else:
            raise ValueError

    mocker.patch("dephon.cli.main.loadfn", side_effect=side_effect)

    parsed_args = parse_args_main(["pc"])
    expected = Namespace(
        ccd=mock_ccd,
        spline_deg=3,
        fig_name="ccd.pdf",
        func=parsed_args.func)
    assert parsed_args == expected

