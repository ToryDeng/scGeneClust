# -*- coding: utf-8 -*-
# @Time : 2022/6/3 17:09
# @Author : Tory Deng
# @File : _utils.py
# @Software: PyCharm
import anndata as ad
import numpy as np
from sklearn.metrics.pairwise import paired_distances
from sklearn.preprocessing import minmax_scale


def compute_gene_closeness(adata: ad.AnnData, centers: np.ndarray) -> np.ndarray:
    """
    Compute the distance of each gene to its cluster mean, min-max normalized in cluster.
    The gene closeness is 1 - normalized distances and named as 'score' in 'adata.var'.

    :param adata: The AnnData object
    :param centers: The cluster means
    :return: distances of all genes
    """
    all_distances = paired_distances(adata.varm['pca'], centers[adata.var['cluster']])
    for gene_cluster in np.unique(adata.var['cluster']):
        gene_cluster_mask = adata.var['cluster'] == gene_cluster
        all_distances[gene_cluster_mask] = 1 - minmax_scale(all_distances[gene_cluster_mask])
    return all_distances
