# -*- coding: utf-8 -*-
# @Time : 2022/6/3 16:27
# @Author : Tory Deng
# @File : preprocessing.py
# @Software: PyCharm
from typing import Literal

import anndata as ad
import scanpy as sc
from loguru import logger


def preprocess(adata: ad.AnnData, version: Literal['fast', 'ps']):
    """
    Preprocess the input AnnData object.

    :param adata: The AnnData object. Must have raw counts in `.X` and have been controlled quality.
    :param version: The version of GeneClust.
    """
    logger.info(f"Preprocessing the data...")

    # normalization
    if version == 'ps':
        # for cell clustering
        adata.layers['X_cell_norm'] = sc.pp.normalize_total(adata, inplace=False)['X']
        adata.layers['X_cell_log'] = sc.pp.log1p(adata.layers['X_cell_norm'], copy=True)  # after log-normalizations
        adata.layers['X_cell_scale'] = sc.pp.scale(adata.layers['X_cell_log'], copy=True)  # after scaling

    # for gene clustering
    adata.layers['X_gene_log'] = sc.pp.log1p(adata.X, copy=True)  # after log-transformation
    adata.layers['X_gene_scale'] = sc.pp.scale(adata.layers['X_gene_log'].T, copy=True).T
    logger.info(f"Data preprocessing finished!")
