# -*- coding: utf-8 -*-
#  Copyright (c) 2023 Kumagai group.
from dataclasses import dataclass

import numpy as np
from monty.json import MSONable
from vise.util.mix_in import ToJsonFileMixIn


@dataclass
class ScreeningCharge(MSONable, ToJsonFileMixIn):
    defect_charge: int
    defect_charge_diff: int
    epsilon_static: np.ndarray
    epsilon_ionic: np.ndarray
    disp: float

    @property
    def epsilon_prime(self):
        ave_static = np.trace(self.epsilon_static) / 3
        ave_ionic = np.trace(self.epsilon_ionic) / 3
        return 1.0 - ave_static / (ave_ionic + ave_static)

    @property
    def charge(self):
        first_term = - self.defect_charge * self.epsilon_prime
        second_term = - self.defect_charge_diff * self.epsilon_prime * self.disp
        return first_term + second_term
