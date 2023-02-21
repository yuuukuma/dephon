# -*- coding: utf-8 -*-
#  Copyright (c) 2023 Kumagai group.

from dephon.util import reduce_wswq


def reduce_wswq_auto(args):
    min_info = args.dephon_init.min_info_from_charge(args.single_ccd.charge)
    band_indices = args.band_indices or min_info.relevant_band_indices
    reduce_wswq(args.wswq, band_indices)