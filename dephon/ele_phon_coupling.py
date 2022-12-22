# -*- coding: utf-8 -*-
#  Copyright (c) 2022 Kumagai group.
from abc import ABC, abstractmethod
from dataclasses import dataclass
from typing import List

from monty.json import MSONable
from nonrad.scaling import sommerfeld_parameter
from vise.util.mix_in import ToJsonFileMixIn

from dephon.enum import Carrier


@dataclass
class CalcWij:
    defect_band_index: int
    band_edge_index: int
    inner_products: List[float]
    Qs: List[float]
    eigenvalue_diff: float

    @property
    def wij(self) -> float:
        # np.polyfit(Q, matels[i, :], 1)[0])
        pass


@dataclass
class EPCoupling:
    charge: int
    capture_type: Carrier
    spin: int
    volume: float
    calc_wif: List[CalcWij]

    def scalings(self, T, m_eff, dielectric_const):
        result = []
        if self.charge:
            result.append(
                SommerfeldScaling(T, self.charge, self.capture_type, m_eff, dielectric_const))
        return result

    # @property
    # def averaged_wif(self):
    #     return float(np.average(self.wif))


# def make_ep_coupling_from_files(charge, )
#     return EPCoupling(charge=)


class Scaling(ABC, MSONable, ToJsonFileMixIn):
    @abstractmethod
    def factor(self):
        pass


@dataclass
class SommerfeldScaling(Scaling):
    T: float
    charge: int
    capture_type: Carrier
    m_eff: float
    dielectric_const: float

    @property
    def factor(self):
        Z = self.charge * self.capture_type.charge
        return sommerfeld_parameter(T=self.T, Z=Z, m_eff=self.m_eff,
                                    eps0=self.dielectric_const)


# @dataclass
# class ChargedSupercellScaling(Scaling):
#     capture_type: Carrier

    # @property
    # def factor(self):
    #     return 0
