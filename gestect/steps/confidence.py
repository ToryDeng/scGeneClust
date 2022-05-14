# -*- coding: utf-8 -*-
# @Time : 2022/5/10 11:40
# @Author : Tory Deng
# @File : confidence.py
# @Software: PyCharm
from queue import Queue
from typing import Literal
from typing import Optional

import igraph
import leidenalg
import numpy as np
import pandas as pd
from scipy.spatial.distance import squareform
from sklearn.metrics import pairwise_distances
from sklearn.metrics.pairwise import pairwise_kernels
from sklearn.mixture import GaussianMixture
from sklearn.neighbors import NearestNeighbors


def find_confident_cells(expr: np.ndarray,
                         n_cell_clusters: int,
                         queue: Queue,
                         how: Literal['single_proba', 'mnn_proba', 'mnn_graph'],
                         k_neignbors: Optional[int] = None,
                         ):
    print('Cell clustering starts...')
    if how == 'single_proba':
        gmm = GaussianMixture(n_components=n_cell_clusters, random_state=2022)
        gmm.fit(expr)
        cell_proba = gmm.predict_proba(expr)
        max_proba = cell_proba.max(1)
        is_confident = max_proba > 0.95
        cluster_label = cell_proba.argmax(1)
    elif how == 'mnn_proba':
        # find confident cells using gmm
        gmm = GaussianMixture(n_components=n_cell_clusters, random_state=2022)
        gmm.fit(expr)
        cell_proba = gmm.predict_proba(expr)
        max_proba = cell_proba.max(1)
        gmm_is_confident = max_proba > 0.95
        cluster_label = cell_proba.argmax(1)
        confident_cluster_label = cluster_label[gmm_is_confident]  # cluster label of confident cells
        gmm_conf_expr = pd.DataFrame(expr[gmm_is_confident])  # expr of confident cells after gmm

        # find mnn in each cluster
        neighbor_k = 10
        is_confident = np.array([False] * len(expr))
        for i in np.unique(cluster_label):
            if pd.value_counts(cluster_label)[i] >= neighbor_k:
                cluster_i = confident_cluster_label == i
                nbrs = NearestNeighbors(n_neighbors=neighbor_k).fit(gmm_conf_expr[cluster_i])
                indices = nbrs.kneighbors(gmm_conf_expr, return_distance=False)
                for j in range(len(indices)):
                    for neighbor in indices[j]:  # find all neighbors of cell j
                        if j in indices[neighbor]:  # if cell j is also the neighbor of its neighbor
                            is_confident[list(gmm_conf_expr[cluster_i].index)[j]] = True
                            is_confident[list(gmm_conf_expr[cluster_i].index)[neighbor]] = True
    elif how == 'mnn_graph':
        cluster_label, is_confident = mnn_graph_clustering(expr, k_neignbors)
    else:
        raise NotImplementedError(f"{how} has not been implemented!")
    queue.put({'cluster': cluster_label, 'confident': is_confident})


def mnn_graph_clustering(cell_embeddings: np.ndarray,
                         k_neighbors: int = 30):
    # TODO: convert similarity matrix to distance matrix
    simi_mtx = pd.DataFrame(data=cell_embeddings.T).corr(method='spearman').values
    ident = np.ones(shape=(simi_mtx.shape[0], 1))
    diag = np.expand_dims(np.diag(simi_mtx), axis=1)
    dis_mtx = np.matmul(ident, diag.T) + np.matmul(diag, ident.T) - 2 * simi_mtx
    # simi_mtx = pairwise_kernels(cell_embeddings, metric='rbf', n_jobs=-1)  # cell similarity matrix
    # create mnn graph
    neigh = NearestNeighbors(n_neighbors=k_neighbors, metric='precomputed', n_jobs=-1)
    neigh.fit(dis_mtx)  # convert similarity matrix to distance matrix
    knn_graph = neigh.kneighbors_graph()
    mnn_dense_graph = [1 if knn_graph[i, j] == 1 and knn_graph[j, i] == 1 else 0
                       for i in range(cell_embeddings.shape[0] - 1) for j in range(i + 1, cell_embeddings.shape[0])]
    mnn_graph = squareform(mnn_dense_graph)
    G = igraph.Graph.Adjacency((mnn_graph > 0).tolist())
    # add edge weights and node labels
    G.es['weight'] = mnn_graph[mnn_graph.nonzero()]

    partition = leidenalg.find_partition(G, partition_type=leidenalg.RBConfigurationVertexPartition, n_iterations=-1,
                                         resolution_parameter=1, seed=2022)
    cluster_label = np.array(partition.membership)
    cluster_count = pd.Series(cluster_label).value_counts(ascending=False)
    conf_cluster_idx = cluster_count[cluster_count > 1].index
    is_confident = np.isin(cluster_label, conf_cluster_idx)
    print('number of cell clusters containing more than 1 cell: ', conf_cluster_idx.shape[0])
    return cluster_label, is_confident
