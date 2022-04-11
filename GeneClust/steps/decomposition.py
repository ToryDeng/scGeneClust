# -*- coding: utf-8 -*-
# @Time : 2022/4/8 21:00
# @Author : Tory Deng
# @File : decomposition.py
# @Software: PyCharm
"""
1. PCA
2. GLM-PCA
3. UMAP
"""
import pandas as pd
import umap
from typing import Literal
import scanpy as sc
import anndata as ad
from glmpca.glmpca import glmpca


def reduce_dimension(
        adata: ad.AnnData,
        dr_method: Literal['pca', 'glm-pca', 'umap'],
        n_componets: int
):
    if dr_method == 'pca':
        embedding = sc.pp.pca(adata.T.X, n_comps=n_componets, copy=True)  # compute PCs on normalized data
    elif dr_method == 'glm-pca':
        embedding = glmpca(adata.raw.X, L=n_componets)["factors"]  # use raw counts
    elif dr_method == 'umap':
        embedding = umap.UMAP(n_components=n_componets).fit_transform(adata.T.X)
    else:
        raise NotImplementedError(f"{dr_method} has not been implemented!")
    return pd.DataFrame(embedding, index=adata.var_names)

