# -*- coding: utf-8 -*-
# @Time : 2022/5/25 12:51
# @Author : Tory Deng
# @File : test_scGeneClust.py
# @Software: PyCharm
from scGeneClust import scGeneClust
from scGeneClust._utils import load_example_adata


def test_scGeneClust():
    adata = load_example_adata(min_cells=2)
    raw_adata = adata.raw.to_adata()
    genes = scGeneClust(raw_adata, mode='fast', n_gene_clusters=200)


