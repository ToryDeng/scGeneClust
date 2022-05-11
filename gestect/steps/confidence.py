# -*- coding: utf-8 -*-
# @Time : 2022/5/10 11:40
# @Author : Tory Deng
# @File : confidence.py
# @Software: PyCharm
import numpy as np
from queue import Queue
from typing import Literal
from sklearn.mixture import GaussianMixture
from sklearn.neighbors import NearestNeighbors
import pandas as pd


def find_confident_cells(expr: np.ndarray,
                         n_cell_clusters: int,
                         queue: Queue,
                         how: Literal['single_proba', 'mnn_proba', 'mnn_graph']
                         ):
    print('Cell clustering starts...')
    if how == 'single_proba':
        gmm = GaussianMixture(n_components=n_cell_clusters)
        gmm.fit(expr)
        cell_proba = gmm.predict_proba(expr)
        max_proba = cell_proba.max(1)
        is_confident = max_proba > 0.95
        cluster_label = cell_proba.argmax(1)
    elif how == 'mnn_proba':
        # find confident cells using gmm
        gmm = GaussianMixture(n_components=n_cell_clusters)
        gmm.fit(expr)
        cell_proba = gmm.predict_proba(expr)
        max_proba = cell_proba.max(1)
        gmm_is_confident = max_proba > 0.95
        cluster_label = cell_proba.argmax(1)
        confident_cluster_label = cluster_label[gmm_is_confident]  # cluster label of confident cells
        gmm_conf_expr = pd.DataFrame(expr)[gmm_is_confident]  # expr of confident cells after gmm

        # find mnn in each cluster

        is_confident = np.array([False] * len(expr))
        for i in np.unique(cluster_label):  # for each cluster
            cluster_i = confident_cluster_label == i
            nbrs = NearestNeighbors(n_neighbors=10).fit(gmm_conf_expr[cluster_i])
            indices = nbrs.kneighbors(gmm_conf_expr, return_distance=False)
            for j in range(len(indices)):
                for neighbor in indices[j]:  # find all neighbor of cell j
                    if j in indices[neighbor]:  # if cell j is also the neighbor of its neighbor
                        is_confident[list(gmm_conf_expr[cluster_i].index)[j]] = True
                        is_confident[list(gmm_conf_expr[cluster_i].index)[neighbor]] = True
        is_confident = np.array(is_confident['is_conf'])
    else:
        raise NotImplementedError(f"{how} has not been implemented!")
    queue.put({'cluster': cluster_label, 'confident': is_confident})
