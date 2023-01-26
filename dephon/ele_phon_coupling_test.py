# -*- coding: utf-8 -*-
#  Copyright (c) 2022 Kumagai group.
from matplotlib import pyplot as plt
from numpy.testing import assert_almost_equal
from vise.tests.helpers.assertion import assert_json_roundtrip

from dephon.ele_phon_coupling import InnerProduct


def test_e_p_coupling_to_json_file(e_p_coupling, tmpdir):
    assert_json_roundtrip(e_p_coupling, tmpdir)


def test_reset_inner_products(e_p_coupling):
    e_p_coupling.reset_inner_prod()
    assert e_p_coupling.e_p_matrix_elements == []


def test_inner_prod_vs_q(e_p_matrix_elem):
    inner_prod_1 = InnerProduct(inner_product=0.2, dQ=0.0)
    inner_prod_2 = InnerProduct(inner_product=1.2, dQ=1.0)
    e_p_matrix_elem.inner_products = [inner_prod_1, inner_prod_2]
    assert e_p_matrix_elem._inner_prod_vs_q == ([0.0, 1.0], [0.2, 1.2])


def test_e_p_matrix_element(e_p_matrix_elem):
    inner_prod_1 = InnerProduct(inner_product=0.2, dQ=0.0)
    inner_prod_2 = InnerProduct(inner_product=1.2, dQ=1.0)
    e_p_matrix_elem.inner_products = [inner_prod_1, inner_prod_2]
    ax = plt.gca()
    assert_almost_equal(e_p_matrix_elem.e_p_matrix_element(ax), 1.0)
    plt.show()


def test_e_p_matrix_element(e_p_coupling):
    print(e_p_coupling)

