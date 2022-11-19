# -*- coding: utf-8 -*-
#  Copyright (c) 2022 Kumagai group.
import shutil
from argparse import Namespace

from dephon.cli.main_function import make_cpp_init


def test_make_cpp_init(test_files, tmpdir):
    print(tmpdir)
    tmpdir.chdir()
    shutil.copytree(test_files / "Na3AgO2" / "Va_O1_0", tmpdir / "Va_O1_0")
    shutil.copytree(test_files / "Na3AgO2" / "Va_O1_1", tmpdir / "Va_O1_1")

    args = Namespace(initial_dir=tmpdir / "Va_O1_0",
                     final_dir=tmpdir / "Va_O1_1",
                     i_to_f_div_ratios=[0.1, 0.2],
                     f_to_i_div_ratios=[0.3, 0.4])
    make_cpp_init(args)

