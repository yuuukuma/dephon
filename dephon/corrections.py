# -*- coding: utf-8 -*-
#  Copyright (c) 2022 Kumagai group.
from dataclasses import dataclass

import yaml
from monty.json import MSONable
from vise.util.mix_in import ToYamlFileMixIn

from dephon.enum import CorrectionType


@dataclass
class DephonCorrection(MSONable, ToYamlFileMixIn):
    energy: float
    correction_type: CorrectionType

    def to_yaml(self):
        return yaml.dump({"energy": self.energy,
                          "correction_type": str(self.correction_type)})

    @classmethod
    def from_yaml(cls, filename: str):
        with open(filename) as file:
            d = yaml.safe_load(file)
        type_ = CorrectionType.from_string(d["correction_type"])
        return cls(energy=d["energy"], correction_type=type_)