# -*- coding: utf-8 -*-
# @Time : 2022/5/25 12:51
# @Author : Tory Deng
# @File : test_GeneClust.py
# @Software: PyCharm
import anndata as ad
from scGeneClust import scGeneClust
from scGeneClust.utils import load_PBMC3k


def test_load_example_adata():
    adata = load_PBMC3k()
    assert adata is not None
    assert isinstance(adata, ad.AnnData)


def test_GeneClust():
    raw_adata = load_PBMC3k()
    genes_fast = scGeneClust(raw_adata, 'fast', 200, random_stat=2022, verbosity=2)
    genes_ps = scGeneClust(raw_adata, 'ps', n_cell_clusters=7, scale=1000, top_percent_relevance=5, random_stat=2022, verbosity=2)
    assert genes_fast is not None and genes_ps is not None


