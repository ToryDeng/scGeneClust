# -*- coding: utf-8 -*-
# @Time : 2022/6/3 16:57
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
        n_comps: int,
        random_stat: Optional[int],
):
    logger.info("Start to reduce dimension...")
    pca_gene = PCA(n_components=n_comps, random_state=random_stat)
    adata.varm['pca'] = pca_gene.fit_transform(adata.layers['X_gene_scale'].T)
    if mode == 'one-way':
        pass
    else:
        pca_cell = PCA(n_components=n_comps, random_state=random_stat)
        adata.obsm['pca'] = pca_cell.fit_transform(adata.layers['X_cell_scale'])
    logger.info("Dimension reduction finished!")