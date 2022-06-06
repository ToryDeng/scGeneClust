# -*- coding: utf-8 -*-
# @Time : 2022/5/25 19:41
# @Author : Tory Deng
# @File : distance.py
# @Software: PyCharm
from typing import Union, Optional

import numpy as np
from scipy.sparse import csr_matrix
from sklearn.metrics.pairwise import pairwise_distances, _VALID_METRICS
from scipy.stats import spearmanr


def compute_pairwise_distances(
        data: np.ndarray,
        metric: str,
        n_jobs: Optional[int] = -1,
        **kwargs
) -> Union[np.ndarray, csr_matrix]:

    if metric in _VALID_METRICS and metric != 'correlation':
        return pairwise_distances(data, metric=metric, n_jobs=n_jobs, **kwargs)
    elif metric == 'spearman':
        corr, _ = spearmanr(data, axis=1)
        np.fill_diagonal(corr, 1)
        return 1 - np.abs(corr)
    elif metric == 'pearson':
        corr = np.corrcoef(data)
        np.fill_diagonal(corr, 1)
        return 1 - np.abs(corr)
    else:
        raise NotImplementedError(f"'{metric}' has not been implemented yet.")







