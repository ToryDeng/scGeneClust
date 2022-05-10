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
from queue import Queue
from sklearn.cluster import AgglomerativeClustering, MeanShift
from sklearn.mixture import GaussianMixture


def clustering_genes(
        adata: ad.AnnData,
        clustering: Literal['agg', 'gmm', 'ms'],
        n_clusters: int
):
    if clustering == 'agg':
        if adata.uns['distance'] == 'euclidean':
            params, to_fit = {'affinity': 'euclidean', 'linkage': 'ward'}, adata.varm[adata.uns['dr_method']]
        else:
            params, to_fit = {'affinity': 'precomputed', 'linkage': 'average'}, adata.varp[adata.uns['distance']]
        adata.var['cluster_label'] = AgglomerativeClustering(n_clusters=n_clusters, **params).fit_predict(to_fit)
    elif clustering == 'gmm':
        model = GaussianMixture(n_components=n_clusters, max_iter=200)
        adata.var['cluster_label'] = model.fit_predict(adata.varm[adata.uns['dr_method']])
        adata.uns['cluster_centroid'] = model.means_
    elif clustering == 'ms':
        model = MeanShift(n_jobs=-1)
        adata.var['cluster_label'] = model.fit_predict(adata.varm[adata.uns['dr_method']])
        adata.uns['cluster_centroid'] = model.cluster_centers_
    else:
        raise NotImplementedError(f"{clustering} has not been implemented!")
    cluster_counts = adata.var['cluster_label'].value_counts(ascending=False)
    print(f"The largest cluster has {cluster_counts.iloc[0]} genes, "
          f"which are {np.round(cluster_counts.iloc[0] / adata.n_vars * 100, decimals=4)}% of all genes.")
