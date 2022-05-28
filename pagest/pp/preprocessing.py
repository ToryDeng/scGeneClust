# -*- coding: utf-8 -*-
# @Time : 2022/5/27 16:43
# @Author : Tory Deng
# @File : preprocessing.py
# @Software: PyCharm
import anndata as ad
import scanpy as sc
import numpy as np
from sklearn.preprocessing import scale
from loguru import logger


def preprocess(adata: ad.AnnData):
    """
    Preprocess the input AnnData object.

    :param adata: The AnnData object. Must have raw counts in `.X` and have been controlled quality.
    :return: None
    """
    assert np.all(adata.X % 1 == 0), ValueError("`adata.X` seems to contain normalized data!")
    logger.info(f"Start to preprocess the data...")
    adata.layers['X_cell'] = adata.X.copy()
    adata.layers['X_gene'] = adata.X.copy()
    # log-normalize cells
    sc.pp.normalize_total(adata, layer='X_cell')
    sc.pp.log1p(adata, layer='X_cell')
    sc.pp.scale(adata, layer='X_cell')
    # scale genes
    sc.pp.log1p(adata, layer='X_gene')
    adata.layers['X_gene'] = sc.pp.scale(adata.layers['X_gene'].T, copy=True).T
