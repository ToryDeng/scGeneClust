# -*- coding: utf-8 -*-
# @Time : 2022/6/3 17:55
# @Author : Tory Deng
# @File : filtering.py
# @Software: PyCharm
import anndata as ad
import numpy as np
from typing import Literal
from loguru import logger


def filter(adata: ad.AnnData, mode: Literal['one-way', 'two-way']) -> ad.AnnData:

    gene_cluster_counts = adata.var['cluster'].value_counts()
    # summary of gene clustering
    logger.debug(f"Total number of gene clusters: {gene_cluster_counts.shape[0]}")
    logger.debug(f"Total number of single gene clusters: {(gene_cluster_counts == 1).sum()}")
    logger.debug(f"Size of top 5 largest gene clusters: \n{gene_cluster_counts.head(5)}")
    logger.debug(f"Size of top 5 smallest gene clusters: \n{gene_cluster_counts.tail(5)}")

    keep_genes = np.isin(adata.var['cluster'], gene_cluster_counts[gene_cluster_counts > 1].index)
    logger.debug(f"Removing {adata.n_vars - keep_genes.sum()} single gene clusters...")

    if mode == 'two-way':
        logger.debug(f"Size of cell clusters: \n{adata.obs['cluster'].value_counts()}")
        keep_cells = adata.obs['highly_confident']
        logger.debug(f"Removing {adata.n_obs - keep_cells.sum()} low-confidence cells...")
    else:
        keep_cells = np.ones(shape=(adata.n_obs,), dtype=bool)

    copied = adata[keep_cells, keep_genes].copy()
    is_constant_genes = np.all(copied.X == copied.X[0, :], axis=0)
    logger.debug(f"Removing {is_constant_genes.sum()} constant genes...")
    copied = copied[:, ~is_constant_genes]

    return copied

