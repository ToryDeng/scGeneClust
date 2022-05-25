# -*- coding: utf-8 -*-
# @Time : 2022/5/21 17:24
# @Author : Tory Deng
# @File : test_decomposition.py
# @Software: PyCharm
import numpy as np

from pagest.tl import reduce_dimension
from pagest.utils import load_example_adata, set_logger


def test_reduce_dimension():
    set_logger(verbose=0)
    adata1 = load_example_adata()
    adata2 = adata1.copy()
    reduce_dimension(adata1, 'one-way', random_stat=42)
    reduce_dimension(adata2, 'two-way', random_stat=42)

    assert 'svd' in adata1.varm and 'svd' not in adata1.obsm
    assert 'svd' in adata2.varm and 'svd' in adata2.obsm
    assert adata2.obsm['svd'].shape[0] == adata2.n_obs and adata2.varm['svd'].shape[0] == adata2.n_vars
    assert np.all(adata1.varm['svd'] == adata2.varm['svd'])



