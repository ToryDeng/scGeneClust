# -*- coding: utf-8 -*-
# @Time : 2022/5/27 16:43
# @Author : Tory Deng
# @File : preprocessing.py
# @Software: PyCharm
import anndata as ad
import scanpy as sc
import numpy as np
from loguru import logger


def preprocess(adata: ad.AnnData):
    """
    Preprocess the input AnnData object.

    :param adata: The AnnData object. Must have raw counts in `.X` and have been controlled quality.
    :return: None
    """
    assert np.all(adata.X % 1 == 0), ValueError("`adata.X` seems to contain normalized data!")

    logger.info(f"Start to preprocess the data...")
    adata.layers['X_cell_ori'] = adata.X.copy()
    adata.layers['X_gene_ori'] = adata.X.copy()
    # for cells
    adata.layers['X_cell_norm'] = sc.pp.normalize_total(adata, layer='X_cell_ori', inplace=False)['X']
    adata.layers['X_cell_log'] = sc.pp.log1p(adata.layers['X_cell_norm'], copy=True)  # after log-normalizations
    adata.layers['X_cell_scale'] = sc.pp.scale(adata.layers['X_cell_log'], copy=True)  # after scaling
    # for genes
    adata.layers['X_gene_log'] = sc.pp.log1p(adata.layers['X_gene_ori'], copy=True)  # after log-transformation
    adata.layers['X_gene_scale'] = sc.pp.scale(adata.layers['X_gene_log'].T, copy=True).T
