# -*- coding: utf-8 -*-
# @Time : 2022/5/21 17:24
# @Author : Tory Deng
# @File : test_decomposition.py
# @Software: PyCharm
import numpy as np

from scGeneClust.pp import reduce_dimension, preprocess
from scGeneClust.utils import load_example_adata, set_logger


def test_reduce_dimension():
    set_logger(verbosity=0)
    adata1 = load_example_adata().raw.to_adata()
    adata2 = adata1.copy()
    preprocess(adata1, 'fast')
    preprocess(adata2, 'ps')
    reduce_dimension(adata1, 'fast', n_comps=50, random_stat=42)
    reduce_dimension(adata2, 'ps', n_comps=50, random_stat=42)

    assert 'pca' in adata1.varm and 'pca' not in adata1.obsm
    assert 'pca' in adata2.varm and 'pca' in adata2.obsm
    assert adata2.obsm['pca'].shape[0] == adata2.n_obs and adata2.varm['pca'].shape[0] == adata2.n_vars
    assert np.all(adata1.varm['pca'] == adata2.varm['pca'])



