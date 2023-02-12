# -*- coding: utf-8 -*-
# @Time : 2022/6/3 16:27
# @Author : Tory Deng
# @File : _preprocessing.py
# @Software: PyCharm
from typing import Literal

import anndata as ad
import numpy as np
import scanpy as sc
from loguru import logger
from scipy.sparse import issparse


def normalize(adata: ad.AnnData, modality: Literal['sc', 'st']):
    """
    Normalize raw counts in `adata.X` using pearson residuals.
    `adata.X` is updated with the normalized counts.

    Parameters
    ----------
    adata : AnnData
        The annotated data matrix of shape `n_obs` × `n_vars`.
        Rows correspond to cells and columns to genes.
    """
    logger.info("Preprocessing data...")
    if modality == 'sc':
        adata.X = adata.X.A if issparse(adata.X) else adata.X
        sc.experimental.pp.normalize_pearson_residuals(adata)
    else:
        sc.pp.normalize_total(adata)
        sc.pp.log1p(adata)


def reduce_dim(adata: ad.AnnData, version: Literal['fast', 'ps'], random_state: int = 0):
    """
    Reduce gene-level dimension (in GeneClust-ps) or cell-level dimension (in GeneClust-fast) of data.
    The cell-level principal components are stored in `adata.varm['X_pca']`.
    For GeneClust-ps, the gene-level principal components are stored in `adata.obsm['X_pca']`.

    Parameters
    ----------
    adata : AnnData
        The annotated data matrix of shape `n_obs` × `n_vars`.
        Rows correspond to cells and columns to genes.
    version : Literal['fast', 'ps']
        The version of GeneClust to be used.
    random_state : int
        Change to use different initial states for the optimization.
    """
    norm_counts = np.asanyarray(adata.X)
    adata.varm['X_pca'] = sc.pp.pca(norm_counts.T, svd_solver='auto', random_state=random_state)
    if version == 'ps':
        adata.obsm['X_pca'] = sc.pp.pca(norm_counts, svd_solver='auto', random_state=random_state)
    logger.info("Data preprocessing done.")
