# -*- coding: utf-8 -*-
#  Copyright (c) 2022 Kumagai group.
from argparse import Namespace
from pathlib import Path

from monty.serialization import loadfn
from pydefect.analyzer.calc_results import CalcResults
from pydefect.analyzer.unitcell import Unitcell
from pymatgen.core import Structure
from vise.input_set.incar import ViseIncar
from vise.input_set.prior_info import PriorInfo

from dephon.cli.main_function import make_ccd_init, make_ccd, plot_ccd, \
    make_ccd_dirs, make_wswq_dirs
from dephon.config_coord import Ccd, CcdInit, ImageStructureInfo, \
    MinimumPointInfo


def test_make_ccd_init(test_files, tmpdir):
    tmpdir.chdir()
    dir_ = test_files / "Na3AgO2"

    args = Namespace(excited_dir=dir_ / "Va_O1_1",
                     ground_dir=dir_ / "Va_O1_0",
                     unitcell=Unitcell.from_yaml(dir_ / "unitcell.yaml"),
                     p_state=loadfn(dir_ / "perfect_band_edge_state.json"))
    make_ccd_init(args)
    actual = loadfn("cc/Va_O1_1to0/ccd_init.json").__str__()
    expected = """name: Va_O1
semiconductor type:  n-type
vbm              1.740  supercell vbm  1.615
cbm              4.705  supercell cbm  4.965
dQ (amu^0.5 Å)   2.242
dR (Å)           0.322
M (amu)         48.367
------------------------------------------------------------
  q   carrier    initial symm    final symm     energy    correction    corrected energy    ZPL
  0     h+e           m              m           6.818         0.000               6.818
  1      e            m              m           5.415         0.239               5.654  1.164
  0                   m              m           3.853         0.000               3.853  1.801"""
    assert actual == expected


def test_make_ccd_dirs(tmpdir, ground_structure, excited_structure,
                       intermediate_structure):
    tmpdir.chdir()
    ccd_init = CcdInit(
        name="test",
        ground_state=MinimumPointInfo(charge=1,
                                      structure=ground_structure,
                                      energy=10.0,
                                      energy_correction=200.0,
                                      carriers=[],
                                      initial_site_symmetry="",
                                      final_site_symmetry="",
                                      parsed_dir=""),
        excited_state=MinimumPointInfo(charge=0,
                                       structure=excited_structure,
                                       energy=10.0,
                                       energy_correction=200.0,
                                       carriers=[],
                                       initial_site_symmetry="",
                                       final_site_symmetry="",
                                       parsed_dir=""),
        vbm=-100.0, cbm=100.0, supercell_vbm=-100.0, supercell_cbm=100.0)

    Path("test").mkdir()
    args = Namespace(ccd_init=ccd_init,
                     e_to_g_div_ratios=[0.5, 1.0],
                     g_to_e_div_ratios=[1.0],
                     calc_dir=Path("test"))
    make_ccd_dirs(args)

    actual = Structure.from_file("excited/disp_0.5/POSCAR")
    assert actual == intermediate_structure

    actual = PriorInfo.load_yaml("excited/disp_0.5/prior_info.yaml")
    expected = PriorInfo(charge=0)
    assert actual == expected

    # dQ = sqrt((0.1*10)**2*6 * Element.H.atomic_mass)
    # dQ / 2 =1.2295974951178128
    actual = loadfn("excited/disp_0.5/image_structure_info.json")
    expected = ImageStructureInfo(dQ=1.2295974951178128, correction=200.0,
                                  correction_type="eFNV")
    assert actual == expected

    actual = loadfn("excited/disp_1.0/image_structure_info.json")
    expected = ImageStructureInfo(dQ=0, correction=200.0,
                                  correction_type="eFNV")
    assert actual == expected

    actual = Structure.from_file("ground/disp_1.0/POSCAR")
    assert actual == excited_structure


# @pytest.fixture
# def ccd():
#     return Ccd(name="test",
#                image_infos=[ImageStructureInfo(-0.2, -1.2, 12.0),
#                                     ImageStructureInfo(0.0, 9.0, 10.0),
#                                     ImageStructureInfo(0.2, -0.8, 8.0)],
#                ground_image_infos=[ImageStructureInfo(-0.2, -2.2, -2.0),
#                                    ImageStructureInfo(0.0, 18.0, 0.0),
#                                    ImageStructureInfo(0.2, -1.8, 2.0)])


def test_make_ccd(tmpdir, mocker, ground_structure):
    print(tmpdir)
    tmpdir.chdir()

    for i in ["ground", "excited"]:
        for disp in [-0.2, 0.2]:
            disp_dir = Path(f"{i}/disp_{disp}")
            disp_dir.mkdir(parents=True)
            image_structure_info = ImageStructureInfo(dQ=10.0 + disp,
                                                      correction=0.1,
                                                      correction_type="eFNV")
            image_structure_info.to_json_file(disp_dir / "image_structure_info.json")
            calc_results = CalcResults(ground_structure, energy=disp,
                                       magnetization=0.0, potentials=[0.0])
            calc_results.to_json_file(str(disp_dir / "calc_results.json"))

    ccd_init = mocker.Mock(spec=CcdInit, autospec=True)
    ccd_init.name = "test"
    ccd_init.vbm = 1.0
    ccd_init.cbm = 2.0
    ccd_init.ground_charge = 1
    ccd_init.excited_charge = 0

    args = Namespace(ccd_init=ccd_init,
                     ground_dirs=[Path("ground/disp_-0.2"), Path("ground/disp_0.2")],
                     excited_dirs=[Path("excited/disp_-0.2"), Path("excited/disp_0.2")],
                     skip_shallow=False)
    make_ccd(args)
    actual: Ccd = loadfn("ccd.json")

    ground_info = [ImageStructureInfo(9.8, -0.2 + 1.0, 0.1, "eFNV"),
                   ImageStructureInfo(10.2, 0.2 + 1.0, 0.1, "eFNV")]
    excited_info = [ImageStructureInfo(9.8, -0.2, 0.1, "eFNV"),
                    ImageStructureInfo(10.2, 0.2, 0.1, "eFNV")]
    ground_info_eg = [ImageStructureInfo(9.8, -0.2 + 2.0, 0.1, "eFNV"),
                      ImageStructureInfo(10.2, 0.2 + 2.0, 0.1, "eFNV")]
    expected = Ccd(name="test",
                   image_infos={"ground": ground_info,
                                "excited + p": excited_info,
                                "ground + p + n": ground_info_eg})
    assert actual == expected


def test_plot_ccd(ccd, tmpdir):
    args = Namespace(ccd=ccd, spline_deg=1, fig_name=tmpdir / "ccd.pdf")
    plot_ccd(args)


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

        incar = ViseIncar({"NSW": 100})
        incar.write_file(Path(f"{state}/disp_-0.2/INCAR"))


    Path(f"excited/disp_-0.2/wswq").mkdir()

    ccd_init = mocker.MagicMock()
    ccd_init.ground_state.dir_path = Path(tmpdir/"ground_original")
    ccd_init.excited_state.dir_path = Path(tmpdir/"excited_original")

    args = Namespace(ground_dirs=[Path(f"ground/disp_-0.2")],
                     excited_dirs=[Path(f"excited/disp_-0.2")],
                     ccd_init=ccd_init)
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
                                    "LWAVE": False,
                                    "LORBIT": "None"})
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