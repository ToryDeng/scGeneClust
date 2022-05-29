# -*- coding: utf-8 -*-
# @Time : 2022/5/21 13:08
# @Author : Tory Deng
# @File : decomposition.py
# @Software: PyCharm
from typing import Literal, Optional

import anndata as ad
from loguru import logger
from sklearn.decomposition import PCA


def reduce_dimension(
        adata: ad.AnnData,
        mode: Literal['one-way', 'two-way'],
        n_comps: int = 50,
        random_stat: Optional[int] = None,
):
    logger.info("Start to reduce dimension...")
    pca_gene = PCA(n_components=n_comps, whiten=True, random_state=random_stat)
    adata.varm['pca'] = pca_gene.fit_transform(adata.layers['X_gene'].T)
    if mode == 'one-way':
        pass
    elif mode == 'two-way':
        pca_cell = PCA(n_components=n_comps, whiten=True, svd_solver='arpack', random_state=random_stat)
        adata.obsm['pca'] = pca_cell.fit_transform(adata.layers['X_cell'])
    else:
        raise ValueError(f"Argument `mode` can only be 'one-way' or 'two-way', not '{mode}'.")
    logger.info("Dimension reduction finished!")
