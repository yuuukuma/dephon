# -*- coding: utf-8 -*-
#  Copyright (c) 2022 Kumagai group.
import pytest

from dephon.config_coord import SinglePointInfo, SingleCcd, Ccd, SingleCcdId
from dephon.dephon_init import MinimumPointInfo, DephonInit
from dephon.enum import Carrier
from dephon.make_config_coord import MakeCcd

band_edges = dict(vbm=1.0, cbm=3.0, supercell_vbm=1.1, supercell_cbm=2.9)


@pytest.fixture
def dephon_init(ground_structure, excited_structure):
    va_o1_0 = MinimumPointInfo(name="Va_O",
                               charge=0,
                               structure=ground_structure,
                               energy=11.0,
                               correction_energy=-1.0,
                               magnetization=0.0,
                               localized_orbitals=[],
                               initial_site_symmetry="2mm",
                               final_site_symmetry="2",
                               parsed_dir="/path/to/Va_O1_0",
                               valence_bands=[[]],
                               conduction_bands=[[]]
                               )
    va_o1_1 = MinimumPointInfo(name="Va_O",
                               charge=1,
                               structure=excited_structure,
                               energy=12.0,
                               correction_energy=-1.0,
                               magnetization=1.0,
                               localized_orbitals=[],
                               initial_site_symmetry="2mm",
                               final_site_symmetry="2",
                               parsed_dir="/path/to/Va_O1_1",
                               valence_bands=[],
                               conduction_bands=[]
                               )
    # transition level = -1.0 from CBM
    return DephonInit(min_points=[va_o1_0, va_o1_1],
                      supercell_volume=10.0,
                      ave_electron_mass=1.0, ave_hole_mass=1.0,
                      ave_static_diele_const=1.0,
                      **band_edges)


common = dict(is_shallow=False, used_for_fitting=True)


def test_make_ccd(excited_structure, ground_structure, dephon_init):
    ground = SingleCcd(SingleCcdId(name="from_0_to_1"),
                       charge=0,
                       point_infos=[SinglePointInfo(dQ=0.0,
                                                    disp_ratio=0.0,
                                                    corrected_energy=-100.0,
                                                    **common),
                                    SinglePointInfo(dQ=10.0,
                                                    disp_ratio=1.0,
                                                    corrected_energy=-90.0,
                                                    **common)])
    excited = SingleCcd(SingleCcdId(name="from_1_to_0"),
                        charge=1,
                        point_infos=[SinglePointInfo(dQ=0.0,
                                                     disp_ratio=0.0,
                                                     corrected_energy=-100.0,
                                                     **common),
                                     SinglePointInfo(dQ=10.0,
                                                     disp_ratio=1.0,
                                                     corrected_energy=-90.0,
                                                     **common)])

    actual = MakeCcd(ground, excited, dephon_init).ccd
    expected = Ccd(name="Va_O", ccds=[
        SingleCcd(SingleCcdId(name="ground"),
                  charge=0,
                  point_infos=[SinglePointInfo(dQ=0.0,
                                               disp_ratio=0.0,
                                               corrected_energy=-100.0,
                                               base_energy=-100.0,
                                               **common),
                               SinglePointInfo(dQ=10.0,
                                               disp_ratio=1.0,
                                               corrected_energy=-90.0,
                                               base_energy=-100.0,
                                               **common)]),
        SingleCcd(SingleCcdId(name="excited", carriers=[Carrier.e]),
                  charge=1,
                  point_infos=[SinglePointInfo(dQ=0.0,
                                               disp_ratio=1.0,
                                               # same as excited
                                               corrected_energy=-87.0,
                                               base_energy=-100.0,
                                               **common),
                               SinglePointInfo(dQ=10.0,
                                               disp_ratio=0.0,
                                               # same as excited
                                               corrected_energy=-97.0,
                                               base_energy=-100.0,
                                               **common)]),
        SingleCcd(SingleCcdId(name="ground", carriers=[Carrier.h,
                                                       Carrier.e]),
                  charge=0,
                  point_infos=[SinglePointInfo(dQ=0.0,
                                               disp_ratio=0.0,
                                               # band gap is added to ground
                                               corrected_energy=-98.0,
                                               base_energy=-100.0,
                                               **common),
                               SinglePointInfo(dQ=10.0,
                                               disp_ratio=1.0,
                                               corrected_energy=-88.0,
                                               base_energy=-100.0,
                                               **common)])])
    assert actual == expected

