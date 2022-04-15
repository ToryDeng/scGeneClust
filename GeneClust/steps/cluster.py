# -*- coding: utf-8 -*-
# @Time : 2022/4/8 21:01
# @Author : Tory Deng
# @File : cluster.py
# @Software: PyCharm
"""
The (hierarchical) clustering
"""

from sklearn.cluster import AgglomerativeClustering
import anndata as ad
import numpy as np


def clustering_genes(adata: ad.AnnData, similarity: str, n_clusters: int):
    adata.var['cluster_labels'] = AgglomerativeClustering(
        affinity='precomputed',
        n_clusters=n_clusters,
        linkage='average'
    ).fit_predict(adata.varp[similarity].T)
    cluster_counts = adata.var['cluster_labels'].value_counts(ascending=False)
    print(f"The largest cluster has {cluster_counts[0]} genes, "
          f"which is {np.round(cluster_counts[0] / adata.n_vars * 100, decimals=4)}% of all genes.")
