# -*- coding: utf-8 -*-
#  Copyright (c) 2022 Kumagai group.
import os
import shutil
from argparse import Namespace
from pathlib import Path

from monty.serialization import loadfn
from pydefect.analyzer.unitcell import Unitcell
from pymatgen.core import Structure
from vise.input_set.incar import ViseIncar
from vise.input_set.prior_info import PriorInfo

from dephon.cli.main_function import make_dephon_init, make_ccd, plot_ccd, \
    make_ccd_dirs, make_wswq_dirs, make_single_point_infos, make_single_ccd, \
    plot_eigenvalues, set_quadratic_fitting_q_range
from dephon.config_coord import SinglePointInfo, SingleCcd, Ccd
from dephon.corrections import DephonCorrection
from dephon.dephon_init import DephonInit, MinimumPointInfo
from dephon.enum import CorrectionType


def test_make_dephon_init(test_files, tmpdir):
    tmpdir.chdir()
    dir_ = test_files / "Na3AgO2"

    args = Namespace(first_dir=dir_ / "Va_O1_1",
                     second_dir=dir_ / "Va_O1_0",
                     unitcell=Unitcell.from_yaml(dir_ / "unitcell.yaml"),
                     p_state=loadfn(dir_ / "perfect_band_edge_state.json"))
    make_dephon_init(args)
    actual = loadfn("cc/Va_O1_1_0/dephon_init.json").__str__()
    expected = """name: Va_O1
vbm              1.740  supercell vbm  1.615
cbm              4.705  supercell cbm  4.965
dQ (amu^0.5 Å)   2.242
dR (Å)           0.322
M (amu)         48.367
------------------------------------------------------------
  q   initial symm    final symm     energy    correction    corrected energy     ZPL
  1        m              m           2.450         0.239               2.689
  0        m              m           3.853         0.000               3.853  -1.164"""
    assert actual == expected


def test_make_ccd_dirs(tmpdir, ground_structure, excited_structure,
                       intermediate_structure):
    print(tmpdir)
    tmpdir.chdir()
    dephon_init = DephonInit(
        defect_name="test",
        states=[MinimumPointInfo(charge=1,
                                 structure=ground_structure,
                                 energy=10.0,
                                 correction_energy=200.0,
                                 initial_site_symmetry="",
                                 final_site_symmetry="",
                                 parsed_dir=""),
                MinimumPointInfo(charge=0,
                                 structure=excited_structure,
                                 energy=10.0,
                                 correction_energy=100.0,
                                 initial_site_symmetry="",
                                 final_site_symmetry="",
                                 parsed_dir="")],
        vbm=-100.0, cbm=100.0, supercell_vbm=-100.0, supercell_cbm=100.0)

    Path("test").mkdir()
    args = Namespace(dephon_init=dephon_init,
                     first_to_second_div_ratios=[0.5, 1.0],
                     second_to_first_div_ratios=[0.0, 1.0],
                     calc_dir=Path("test"))
    make_ccd_dirs(args)

    actual = Structure.from_file("from_1_to_0/disp_0.5/POSCAR")
    assert actual == intermediate_structure

    actual = PriorInfo.load_yaml("from_1_to_0/disp_0.5/prior_info.yaml")
    expected = PriorInfo(charge=1)
    assert actual == expected

    # dQ = sqrt((0.1*10)**2*6 * Element.H.atomic_mass)
    # dQ / 2 =1.2295974951178128
    actual = loadfn("from_1_to_0/disp_0.5/single_point_info.json")
    expected = SinglePointInfo(dQ=1.2295974951178128, disp_ratio=0.5)
    assert actual == expected

    actual = loadfn("from_0_to_1/disp_0.0/single_point_info.json")
    expected = SinglePointInfo(dQ=0.0, disp_ratio=0.0)
    assert actual == expected

    actual = Structure.from_file("from_1_to_0/disp_1.0/POSCAR")
    assert actual == excited_structure

    actual = DephonCorrection.from_yaml("from_1_to_0/disp_1.0/dephon_correction.yaml")
    expected = DephonCorrection(200.0, CorrectionType.extended_FNV)
    assert actual == expected


# @pytest.fixture
# def ccd():
#     return Ccd(name="test",
#                image_infos=[SinglePointInfo(-0.2, -1.2, 12.0),
#                                     SinglePointInfo(0.0, 9.0, 10.0),
#                                     SinglePointInfo(0.2, -0.8, 8.0)],
#                ground_image_infos=[SinglePointInfo(-0.2, -2.2, -2.0),
#                                    SinglePointInfo(0.0, 18.0, 0.0),
#                                    SinglePointInfo(0.2, -1.8, 2.0)])
#
#


def test_make_single_point_infos(test_files, tmpdir):
    tmpdir.chdir()
    src = Path(test_files / "NaP/Va_P1_-1_0/from_0_to_-1_before_make_single_point_infos/disp_0.0")
    print(src)
    shutil.copytree(src, Path.cwd() / "disp_0.0")
    args = Namespace(dirs=[Path("disp_0.0")])
    make_single_point_infos(args)

    actual = loadfn("disp_0.0/single_point_info.json")
    expected = SinglePointInfo(dQ=0.0,
                               disp_ratio=0.0,
                               corrected_energy=-2223.75521961,
                               is_shallow=False,
                               correction_method=CorrectionType.extended_FNV)
    assert actual == expected


def test_make_single_ccd(test_files, tmpdir):
    tmpdir.chdir()
    print(tmpdir)
    src = Path(test_files / "NaP/Va_P1_-1_0/from_0_to_-1_after_make_single_point_infos")
    shutil.copytree(src, Path.cwd() / "from_0_to_-1")
    os.chdir(Path("from_0_to_-1"))
    args = Namespace(dirs=[Path("disp_0.0")])
    make_single_ccd(args)

    actual = loadfn("single_ccd.json")

    point_info_disp = \
        SinglePointInfo(dQ=0.0,
                        disp_ratio=0.0,
                        corrected_energy=-2223.75521961,
                        is_shallow=False,
                        correction_method=CorrectionType.extended_FNV)
    expected = SingleCcd(name="from_0_to_-1",
                         charge=0,
                         point_infos=[point_info_disp])
    assert actual == expected


def test_make_ccd(test_files, tmpdir):
    tmpdir.chdir()
    va_p1 = Path(test_files) / "NaP/Va_P1_-1_0"
    ground_ccd = loadfn(va_p1 / "from_-1_to_0_after_make_single_point_infos/single_ccd.json")
    excited_ccd = loadfn(va_p1 / "from_0_to_-1_after_make_single_point_infos/single_ccd.json")
    dephon_init = loadfn(va_p1 / "dephon_init.json")
    args = Namespace(ground_ccd=ground_ccd, excited_ccd=excited_ccd,
                     dephon_init=dephon_init)
    make_ccd(args)


def test_set_quadratic_fitting_q_range(ccd, tmpdir):
    tmpdir.chdir()
    args = Namespace(ccd=ccd, single_ccd_name="ground", q_range=[-1.0, 1.0])
    set_quadratic_fitting_q_range(args)
    ccd: Ccd = loadfn("ccd.json")
    assert ccd.ccds[1].point_infos[1].used_for_fitting is True


def test_plot_ccd(ccd, tmpdir):
    args = Namespace(ccd=ccd, fig_name=tmpdir / "ccd.pdf")
    plot_ccd(args)


def test_plot_eigenvalues(test_files, tmpdir):
    tmpdir.chdir()
    va_p1 = Path(test_files) / "NaP/Va_P1_-1_0"
    dephon_init = loadfn(va_p1 / "dephon_init.json")
    dir_ = va_p1 / "from_0_to_-1_after_make_single_point_infos"
    args = Namespace(dirs=[dir_ / "disp_0.0", dir_ / "disp_1.0"],
                     dephon_init=dephon_init)
    plot_eigenvalues(args)


def test_make_wswq_dirs(tmpdir, mocker):

    print(tmpdir)
    tmpdir.chdir()

    for state in ["ground", "excited"]:
        Path(f"{state}_original").mkdir(parents=True)
        Path(f"{state}/disp_-0.2").mkdir(parents=True)
        Path(f"{state}_original/WAVECAR").write_text("wave")
        Path(f"{state}/disp_-0.2/KPOINTS").write_text("kpoints")
        Path(f"{state}/disp_-0.2/POSCAR").write_text("poscar")
        Path(f"{state}/disp_-0.2/POTCAR").write_text("potcar")
        Path(f"{state}/disp_-0.2/WAVECAR").write_text("qqq")

        incar = ViseIncar({"NSW": 100, "LORBIT": 11})
        incar.write_file(Path(f"{state}/disp_-0.2/INCAR"))

    Path(f"excited/disp_-0.2/wswq").mkdir()

    dephon_init = mocker.MagicMock()
    dephon_init.single_ccd.dir_path = Path(tmpdir / "ground_original")
    dephon_init.excited_state.dir_path = Path(tmpdir/"excited_original")

    args = Namespace(ground_dirs=[Path(f"ground/disp_-0.2")],
                     excited_dirs=[Path(f"excited/disp_-0.2")],
                     dephon_init=dephon_init)
    make_wswq_dirs(args)

    for state in ["ground"]:
        wswq_dir = Path(f"{state}/disp_-0.2/wswq/")
        assert Path(wswq_dir/"KPOINTS").read_text() == "kpoints"
        assert Path(wswq_dir/"POSCAR").read_text() == "poscar"
        assert Path(wswq_dir/"POTCAR").read_text() == "potcar"
        assert Path(wswq_dir/"WAVECAR.qqq").read_text() == "qqq"
        assert Path(wswq_dir/"WAVECAR").read_text() == "wave"

        actual_incar = ViseIncar.from_file(wswq_dir/"INCAR")
        expected_incar = ViseIncar({"NSW": 100,
                                    "ALGO": "None",
                                    "LWSWQ": True,
                                    "NELM": 1,
                                    "LWAVE": False})
        assert actual_incar == expected_incar

        assert Path(wswq_dir/"KPOINTS").is_symlink()
        assert Path(wswq_dir/"POSCAR").is_symlink()
        assert Path(wswq_dir/"POTCAR").is_symlink()
        assert Path(wswq_dir/"WAVECAR.qqq").is_symlink()
        assert Path(wswq_dir/"WAVECAR").is_symlink()

    assert Path("excited/disp_-0.2/wswq/KPOINTS").exists() is False


"""
TODO:
. Check if electronic SCF are converged.
. Check if displace_ratio=1 structure is the same as the counterpart.
  Consider the energy correction
"""