# -*- coding: utf-8 -*-
# @Time : 2022/6/3 16:27
# @Author : Tory Deng
# @File : preprocessing.py
# @Software: PyCharm
from typing import Literal

import anndata as ad
import scanpy as sc
from loguru import logger


def preprocess(adata: ad.AnnData, mode: Literal['fast', 'hc']):
    """
    Preprocess the input AnnData object.

    :param adata: The AnnData object. Must have raw counts in `.X` and have been controlled quality
    :param mode:
    :return: None
    """
    logger.info(f"Start to preprocess the data...")

    # normalization
    if mode == 'hc':
        # for cells
        adata.layers['X_cell_norm'] = sc.pp.normalize_total(adata, inplace=False)['X']
        adata.layers['X_cell_log'] = sc.pp.log1p(adata.layers['X_cell_norm'], copy=True)  # after log-normalizations
        adata.layers['X_cell_scale'] = sc.pp.scale(adata.layers['X_cell_log'], copy=True)  # after scaling

    # for genes
    adata.layers['X_gene_log'] = sc.pp.log1p(adata.X, copy=True)  # after log-transformation
    adata.layers['X_gene_scale'] = sc.pp.scale(adata.layers['X_gene_log'].T, copy=True).T
    logger.info(f"Data preprocessing finished!")

