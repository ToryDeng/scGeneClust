# -*- coding: utf-8 -*-
# @Time : 2022/5/10 11:49
# @Author : Tory Deng
# @File : gene_clustering.py
# @Software: PyCharm
from queue import Queue

import numpy as np
from scipy.cluster.hierarchy import linkage, cophenet, fcluster
from scipy.spatial.distance import pdist
from sklearn.cluster import MeanShift
from sklearn.mixture import GaussianMixture


def cluster_genes(expr: np.ndarray,
                  method: str,
                  n_gene_clusters: int = None,
                  queue: Queue = None):
    print('Gene clustering starts...')
    if method == 'ms':
        ms = MeanShift()
        cluster_labels = ms.fit_predict(expr)
    elif method == 'agg':
        hierarchy = linkage(expr, method='ward')
        c, coph_dists = cophenet(hierarchy, pdist(expr))
        print("coph: ", c)
        cluster_labels = fcluster(hierarchy, n_gene_clusters, criterion='maxclust')
    elif method == 'gmm':
        gmm = GaussianMixture(n_components=n_gene_clusters, random_state=2022)
        cluster_labels = gmm.fit_predict(expr)
    else:
        raise NotImplementedError(f"'{method}' has not been implemented yet!")
    queue.put(cluster_labels)
