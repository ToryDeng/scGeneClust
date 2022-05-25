# -*- coding: utf-8 -*-
# @Time : 2022/5/23 12:13
# @Author : Tory Deng
# @File : filter.py
# @Software: PyCharm
from typing import Literal

import anndata as ad
import numpy as np
from loguru import logger


def filter_adata(adata: ad.AnnData, mode: Literal['one-way', 'two-way'], quantile: float = 0.1) -> ad.AnnData:
    """
    Filter out low-quality gene clusters (and low-confidence cells if `mode`='two-way') on copied AnnData object.

    :param adata: The AnnData object
    :param mode: `one-way` only considers patterns in genes; `two-way` considers patterns both in cells and genes
    :param quantile: The quantile of number of gene clusters to compute
    :return: The copied and filtered AnnData object
    """
    cluster_score = adata.var.loc[:, ('cluster', 'cluster_score')].drop_duplicates()['cluster_score']
    bound = cluster_score.quantile(q=quantile)
    keep_genes = adata.var['cluster_score'] > bound
    logger.debug(f"Removing {np.sum(cluster_score < bound)} gene clusters ({adata.n_vars - keep_genes.sum()} genes)")

    if mode == 'two-way':
        keep_cells = adata.obs['highly_confident']
        logger.debug(f"Removing {adata.n_obs - keep_cells.sum()} low-confidence cells")
    else:
        keep_cells = np.ones(shape=(adata.n_obs, ), dtype=bool)
    copied = adata[keep_cells, keep_genes].copy()
    is_constant_genes = np.all(copied.X == copied.X[0, :], axis=0)
    logger.debug(f"Removing {is_constant_genes.sum()} constant genes...")
    copied = copied[:, ~is_constant_genes]

    return copied







