# -*- coding: utf-8 -*-
# @Time : 2022/5/21 13:08
# @Author : Tory Deng
# @File : decomposition.py
# @Software: PyCharm
from typing import Literal, Optional

import anndata as ad
from loguru import logger
from sklearn.decomposition import TruncatedSVD


def reduce_dimension(
        adata: ad.AnnData,
        mode: Literal['one-way', 'two-way'],
        layer: Optional[str] = None,
        n_comps: int = 50,
        random_stat: Optional[int] = None,
):
    logger.info("Start to reduce dimension...")
    expression_matrix = adata.layers[layer] if layer is not None else adata.X
    decomp_model = TruncatedSVD(n_components=n_comps, n_iter=10, random_state=random_stat)
    decomp_model.fit(expression_matrix)
    if mode == 'one-way':
        adata.varm['svd'] = decomp_model.components_.T
    elif mode == 'two-way':
        adata.varm['svd'], adata.obsm['svd'] = decomp_model.components_.T, decomp_model.transform(expression_matrix)
    else:
        raise ValueError(f"Argument `mode` can only be 'one-way' or 'two-way', not '{mode}'.")
    logger.info("Dimension reduction finished!")
