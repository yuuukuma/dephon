# -*- coding: utf-8 -*-
#  Copyright (c) 2022 Kumagai group.
import re
from pathlib import Path
from typing import List, Dict, Tuple

from pymatgen.electronic_structure.core import Spin


def spin_to_idx(spin: Spin, count_from_1=False) -> int:
    result = 0 if spin == Spin.up else 1
    if count_from_1:
        result += 1
    return result


def idx_to_spin(idx: int):
    return Spin.up if idx == 0 else Spin.down


def reduce_wswq(filename: Path, orbitals: Dict[Tuple[int, int], List[int]]):
    lines = []
    with filename.open(mode='r') as f:
        for line in f:
            header = re.search(r'\s*spin=(\d+), kpoint=\s*(\d+)', line)
            if header:
                spin_kpt = (int(header.group(1)), int(header.group(2)))
                if spin_kpt in orbitals:
                    lines.append(line)
                    orbs = orbitals[spin_kpt]
            else:
                _append_ij(line, orbs, lines)

    filename.write_text("".join(lines))


def _append_ij(line, orbitals, out):
    data = re.search(r'i=\s*(\d+), '
                     r'' r'j=\s*(\d+)\s*:\s*([0-9\-.]+)\s+([0-9\-.]+)',
                     line)
    i, j = int(data.group(1)), int(data.group(2))
    if i in orbitals and j in orbitals:
        out.append(line)
