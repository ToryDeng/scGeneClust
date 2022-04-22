# -*- coding: utf-8 -*-
# @Time : 2022/4/8 21:01
# @Author : Tory Deng
# @File : cluster.py
# @Software: PyCharm
"""
The (hierarchical) clustering
"""

from typing import Literal

import anndata as ad
import numpy as np
from sklearn.cluster import AgglomerativeClustering
from sklearn.mixture import GaussianMixture


def clustering_genes(
        adata: ad.AnnData,
        clustering: Literal['agglomerative', 'gmm'],
        n_clusters: int
):
    if clustering == 'agglomerative':
        adata.var['cluster_label'] = AgglomerativeClustering(
            affinity='precomputed',
            n_clusters=n_clusters,
            linkage='average'
        ).fit_predict(adata.varp[adata.uns['distance']])
    elif clustering == 'gmm':
        adata.var['cluster_label'] = GaussianMixture(
            n_components=n_clusters,
            max_iter=200
        ).fit_predict(adata.varm[adata.uns['dr_method']])
    else:
        raise NotImplementedError(f"{clustering} has not been implemented!")
    cluster_counts = adata.var['cluster_label'].value_counts(ascending=False)
    print(f"The largest cluster has {cluster_counts.iloc[0]} genes, "
          f"which are {np.round(cluster_counts.iloc[0] / adata.n_vars * 100, decimals=4)}% of all genes.")
