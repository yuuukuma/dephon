# -*- coding: utf-8 -*-
#  Copyright (c) 2022 Kumagai group.
from argparse import Namespace
from pathlib import Path

from pydefect.analyzer.band_edge_states import PerfectBandEdgeState

from dephon.cli.main import parse_args_main
from dephon.config_coord import Ccd
from dephon.dephon_init import DephonInit


def loadfn_effect(d: dict):
    def side_effect(filename):
        try:
            return d[filename]
        except KeyError:
            raise ValueError
    return side_effect


def test_main_make_dephon_init(mocker):
    mock_unitcell = mocker.patch("pydefect.cli.main.Unitcell")
    mock_p_band_edge_state = mocker.Mock(spec=PerfectBandEdgeState, autospec=True)
    side_effect = loadfn_effect(
        {"perfect_band_edge_state.json": mock_p_band_edge_state})
    mocker.patch("pydefect.cli.main.loadfn", side_effect=side_effect)
    parsed_args = parse_args_main(["mdi",
                                   "-fd", "Va_O1_0",
                                   "-sd", "Va_O1_1",
                                   "-u", "unitcell.yaml",
                                   "-pbes", "perfect_band_edge_state.json"
                                   ])
    expected = Namespace(
        first_dir=Path("Va_O1_0"),
        second_dir=Path("Va_O1_1"),
        unitcell=mock_unitcell.from_yaml.return_value,
        p_state=mock_p_band_edge_state,
        func=parsed_args.func)
    assert parsed_args == expected
    mock_unitcell.from_yaml.assert_called_once_with("unitcell.yaml")


def test_main_make_dirs_wo_args(mocker):
    mock_dephon_init = mocker.Mock(spec=DephonInit, autospec=True)
    side_effect = loadfn_effect({"dephon_init.json": mock_dephon_init})
    mocker.patch("dephon.cli.main.loadfn", side_effect=side_effect)

    parsed_args = parse_args_main(["mcd"])
    expected = Namespace(
        dephon_init=mock_dephon_init,
        first_to_second_div_ratios=[-0.4, -0.2, 0.0, 0.2, 0.4, 0.6, 0.8],
        second_to_first_div_ratios=[-0.4, -0.2, 0.0, 0.2, 0.4, 0.6, 0.8],
        calc_dir=Path.cwd(),
        func=parsed_args.func)
    assert parsed_args == expected


def test_main_make_ccd_wo_args(mocker):
    mock_dephon_init = mocker.Mock(spec=DephonInit, autospec=True)
    side_effect = loadfn_effect({"dephon_init.json": mock_dephon_init})

    mocker.patch("dephon.cli.main.loadfn", side_effect=side_effect)

    parsed_args = parse_args_main(["mc", "-g", "Va_O1_4", "-e", "Va_O1_3"])
    expected = Namespace(
        dephon_init=mock_dephon_init,
        ground_dirs=[Path("Va_O1_4")],
        excited_dirs=[Path("Va_O1_3")],
        skip_shallow=False,
        func=parsed_args.func)
    assert parsed_args == expected


def test_main_plot_ccd_wo_args(mocker):
    mock_ccd = mocker.Mock(spec=Ccd, autospec=True)
    side_effect = loadfn_effect({"ccd.json": mock_ccd})
    mocker.patch("dephon.cli.main.loadfn", side_effect=side_effect)

    parsed_args = parse_args_main(["pc"])
    expected = Namespace(
        ccd=mock_ccd,
        fig_name="ccd.pdf",
        func=parsed_args.func)
    assert parsed_args == expected

