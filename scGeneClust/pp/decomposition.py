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
        version: Literal['fast', 'ps'],
        random_stat: Optional[int],
):
    """

    :param adata: The annotated data matrix.
    :param version: Version of GeneClust.
    :param random_stat: Change to use different initial states for the optimization.
    """
    logger.info("Reducing data dimension...")
    if version == 'fast':
        pca_gene = PCA(n_components=50, random_state=random_stat)
        adata.varm['pca'] = pca_gene.fit_transform(adata.layers['X_gene_scale'].T)
    else:
        pca_cell = PCA(n_components=50, random_state=random_stat)
        adata.obsm['pca'] = pca_cell.fit_transform(adata.layers['X_cell_scale'])
    logger.info("Dimension reduction finished!")
