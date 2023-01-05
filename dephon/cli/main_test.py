# -*- coding: utf-8 -*-
#  Copyright (c) 2022 Kumagai group.
from argparse import Namespace
from pathlib import Path

from pydefect.analyzer.band_edge_states import PerfectBandEdgeState

from dephon.cli.main import parse_args_main
from dephon.cli.main_function import make_single_point_infos, make_ccd, \
    plot_eigenvalues, make_single_ccd
from dephon.config_coord import Ccd, SingleCcd
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


def test_main_make_dirs(mocker):
    mock_dephon_init = mocker.Mock(spec=DephonInit, autospec=True)
    side_effect = loadfn_effect({"dephon_init.json": mock_dephon_init})
    mocker.patch("dephon.cli.main.loadfn", side_effect=side_effect)

    parsed_args = parse_args_main(["mcd"])
    expected = Namespace(
        dephon_init=mock_dephon_init,
        first_to_second_div_ratios=[-0.2, 0.0, 0.2, 0.4, 0.6, 0.8],
        second_to_first_div_ratios=[-0.2, 0.0, 0.2, 0.4, 0.6, 0.8],
        calc_dir=Path.cwd(),
        func=parsed_args.func)
    assert parsed_args == expected

    parsed_args = parse_args_main(["mcd", "-fsr", "0.0", "-sfr", "0.1",
                                   "-d", "dirname"])
    expected = Namespace(
        dephon_init=mock_dephon_init,
        first_to_second_div_ratios=[0.0],
        second_to_first_div_ratios=[0.1],
        calc_dir=Path("dirname"),
        func=parsed_args.func)
    assert parsed_args == expected


def test_main_make_single_point_infos():
    parsed_args = parse_args_main(["mspi", "-d", "disp_0.0"])
    expected = Namespace(dirs=[Path("disp_0.0")], func=make_single_point_infos)
    assert parsed_args == expected


def test_main_make_single_ccd():
    parsed_args = parse_args_main(["msc", "-d", "disp_0.0"])
    expected = Namespace(dirs=[Path("disp_0.0")], func=make_single_ccd)
    assert parsed_args == expected


def test_main_make_ccd_wo_args(mocker):
    mock_dephon_init = mocker.Mock(spec=DephonInit, autospec=True)
    mock_ground_single_ccd = mocker.Mock(spec=SingleCcd, autospec=True)
    mock_excited_single_ccd = mocker.Mock(spec=SingleCcd, autospec=True)
    side_effect = loadfn_effect(
        {"dephon_init.json": mock_dephon_init,
         "ground_single_ccd.json": mock_ground_single_ccd,
         "excited_single_ccd.json": mock_excited_single_ccd})

    mocker.patch("dephon.cli.main.loadfn", side_effect=side_effect)

    parsed_args = parse_args_main(["mc",
                                   "-g", "ground_single_ccd.json",
                                   "-e", "excited_single_ccd.json"])
    expected = Namespace(
        dephon_init=mock_dephon_init,
        ground_ccd=mock_ground_single_ccd,
        excited_ccd=mock_excited_single_ccd,
        func=make_ccd)
    assert parsed_args == expected


def test_main_set_quadratic_fitting_q_range(mocker):
    mock_ccd = mocker.Mock(spec=Ccd, autospec=True)
    side_effect = loadfn_effect({"ccd.json": mock_ccd})
    mocker.patch("dephon.cli.main.loadfn", side_effect=side_effect)

    parsed_args = parse_args_main(
        ["sfr", "--single_ccd_name", "ground", "--q_range", "-1.0", "1.0"])
    expected = Namespace(
        ccd=mock_ccd,
        single_ccd_name="ground",
        q_range=[-1.0, 1.0],
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
        q_range=None,
        func=parsed_args.func)
    assert parsed_args == expected


def test_main_plot_eigenvalues(mocker):
    mock_dephon_init = mocker.Mock(spec=DephonInit, autospec=True)
    side_effect = loadfn_effect({"dephon_init.json": mock_dephon_init})
    mocker.patch("dephon.cli.main.loadfn", side_effect=side_effect)

    parsed_args = parse_args_main(["pe", "-d", "disp_0.0"])
    expected = Namespace(
        dirs=[Path("disp_0.0")],
        dephon_init=mock_dephon_init,
        func=plot_eigenvalues)
    assert parsed_args == expected

