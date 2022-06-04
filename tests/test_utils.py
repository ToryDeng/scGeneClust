# -*- coding: utf-8 -*-
# @Time : 2022/5/21 17:27
# @Author : Tory Deng
# @File : test_utils.py
# @Software: PyCharm
import anndata as ad

from pagest.utils import load_example_adata


def test_load_example_adata():
    adata = load_example_adata()
    assert isinstance(adata, ad.AnnData)
    assert adata.raw is not None
