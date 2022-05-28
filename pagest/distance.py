# -*- coding: utf-8 -*-
# @Time : 2022/5/25 19:41
# @Author : Tory Deng
# @File : distance.py
# @Software: PyCharm
from typing import Union, Optional

import numpy as np
from scipy.sparse import csr_matrix
from sklearn.metrics.pairwise import pairwise_distances, _VALID_METRICS
from scipy.spatial.distance import squareform


def compute_pairwise_distances(
        X: Union[csr_matrix, np.ndarray],
        metric: str,
        n_jobs: Optional[int] = None,
        **kwargs
):
    """
    Compute pairwise distances between cells or genes.

    :param X: The feature array, where rows represent cells (genes), and columns represent components
    :param metric: The distance metric
    :param n_jobs: The number of jobs to use for the computation
    :param kwargs: Other keyword arguments
    :return:
    """
    if metric in _VALID_METRICS:
        return pairwise_distances(X=X, metric=metric, n_jobs=n_jobs, **kwargs)
    elif metric == 'rho_p':
        data = X.T + 1 / (X.shape[0] ** 2)
        var_log = np.var(np.log(data), axis=0)
        var_mtx = [var_log[i] + var_log[j] for i in range(data.shape[1] - 1) for j in range(i + 1, (data.shape[1]))]
        var_mtx = squareform(var_mtx) + np.diag(2 * var_log)
        return 2 * np.cov(np.log(data.T), ddof=0) / var_mtx
    elif metric == 'phi_s':
        data = X.T + 1 / X.shape[0] ** 2
        var_log = np.var(np.log(data), axis=0)
        phi_s = [np.var(np.log(data[:, i]) - np.log(data[:, j])) / var_log[i] for i in range(data.shape[1]) for j in
                 range(data.shape[1])]
        phi_s = np.array(phi_s).reshape(data.shape[1], data.shape[1])
        return phi_s
    else:
        raise NotImplementedError(f"`{metric}` has not been implemented yet.")



