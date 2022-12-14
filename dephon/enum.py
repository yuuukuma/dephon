# -*- coding: utf-8 -*-
#  Copyright (c) 2022 Kumagai group.
from monty.json import MSONable
from vise.util.enum import ExtendedEnum


class Carrier(MSONable, ExtendedEnum):
    electron = "e"
    hole = "h"

    @property
    def charge(self):
        if self == self.electron:
            return -1
        elif self == self.hole:
            return 1
        else:
            raise ValueError

    @classmethod
    def from_carrier_charge(cls, carrier_charge):
        if carrier_charge == 1:
            return cls.hole
        elif carrier_charge == -1:
            return cls.electron
        raise ValueError


class CorrectionType(MSONable, ExtendedEnum):
    extended_FNV = "extended FNV"
    no_correction = "no correction"


class BandEdge(MSONable, ExtendedEnum):
    vbm = "vbm"
    cbm = "cbm"

