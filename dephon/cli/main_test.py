# -*- coding: utf-8 -*-
#  Copyright (c) 2022 Kumagai group.
from argparse import Namespace
from pathlib import Path

from pydefect.analyzer.band_edge_states import PerfectBandEdgeState

from dephon.cli.main import parse_args_main
from dephon.config_coord import CcdInit, Ccd


def test_main_make_ccd_init(mocker):
    mock_unitcell = mocker.patch("pydefect.cli.main.Unitcell")
    mock_p_band_edge_state = mocker.Mock(spec=PerfectBandEdgeState, autospec=True)

    def side_effect(filename):
        if filename == "perfect_band_edge_state.json":
            return mock_p_band_edge_state
        else:
            raise ValueError

    mocker.patch("pydefect.cli.main.loadfn", side_effect=side_effect)

    parsed_args = parse_args_main(["mci",
                                   "-ed", "Va_O1_0",
                                   "-gd", "Va_O1_1",
                                   "-u", "unitcell.yaml",
                                   "-pbes", "perfect_band_edge_state.json"
                                   ])
    expected = Namespace(
        excited_dir=Path("Va_O1_0"),
        ground_dir=Path("Va_O1_1"),
        unitcell=mock_unitcell.from_yaml.return_value,
        p_state=mock_p_band_edge_state,
        func=parsed_args.func)
    assert parsed_args == expected
    mock_unitcell.from_yaml.assert_called_once_with("unitcell.yaml")


def test_main_make_dirs_wo_args(mocker):
    mock_ccd_init = mocker.Mock(spec=CcdInit, autospec=True)

    def side_effect(filename):
        if filename == "ccd_init.json":
            return mock_ccd_init
        else:
            raise ValueError

    mocker.patch("dephon.cli.main.loadfn", side_effect=side_effect)

    parsed_args = parse_args_main(["mcd"])
    expected = Namespace(
        ccd_init=mock_ccd_init,
        e_to_g_div_ratios=[-0.4, -0.2, 0.2, 0.4, 0.6, 0.8],
        g_to_e_div_ratios=[-0.4, -0.2, 0.2, 0.4, 0.6, 0.8],
        calc_dir=Path.cwd(),
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

    parsed_args = parse_args_main(["mc", "-g", "Va_O1_4", "-e", "Va_O1_3"])
    expected = Namespace(
        ccd_init=mock_ccd_init,
        ground_dirs=[Path("Va_O1_4")],
        excited_dirs=[Path("Va_O1_3")],
        skip_shallow=False,
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

