# -*- coding: utf-8 -*-
# @Time : 2022/5/10 11:40
# @Author : Tory Deng
# @File : confidence.py
# @Software: PyCharm
import numpy as np
from queue import Queue
from typing import Literal
from sklearn.mixture import GaussianMixture


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
    else:
        raise NotImplementedError(f"{how} has not been implemented!")
    queue.put({'cluster': cluster_label, 'confident': is_confident})