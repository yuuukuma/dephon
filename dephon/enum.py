# -*- coding: utf-8 -*-
#  Copyright (c) 2022 Kumagai group.
from monty.json import MSONable
from vise.util.enum import ExtendedEnum


class Carrier(MSONable, ExtendedEnum):
    electron = "e"
    hole = "h"

    @property
    def charge(self):
        if self.electron:
            return -1
        else:
            return 1

    @classmethod
    def from_carrier_charge(cls, carrier_charge):
        if carrier_charge == 1:
            return cls.hole
        elif carrier_charge == -1:
            return cls.electron
        raise ValueError


class BandEdge(MSONable, ExtendedEnum):
    vbm = "vbm"
    cbm = "cbm"

