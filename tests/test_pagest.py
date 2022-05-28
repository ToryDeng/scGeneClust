# -*- coding: utf-8 -*-
# @Time : 2022/5/25 12:51
# @Author : Tory Deng
# @File : test_pagest.py
# @Software: PyCharm
from pagest import pagest
from pagest.utils import load_example_adata


def test_pagest():
    adata = load_example_adata(min_cells=2)
    genes = pagest(adata.raw.to_adata(), n_features=500, mode='one-way', gene_distance='euclidean', verbose=2)
    genes = pagest(adata.raw.to_adata(), n_features=500, mode='two-way', gene_distance='euclidean', verbose=2)
