# -*- coding: utf-8 -*-
# @Time : 2022/5/25 19:41
# @Author : Tory Deng
# @File : distance.py
# @Software: PyCharm
from typing import Union, Optional

import numpy as np
from scipy.sparse import csr_matrix
from sklearn.metrics.pairwise import pairwise_distances, _VALID_METRICS


def compute_pairwise_distances(
        X: Union[csr_matrix, np.ndarray],
        metric: str,
        n_jobs: Optional[int] = None,
        **kwargs
):
    if metric in _VALID_METRICS:
        return pairwise_distances(X=X, metric=metric, n_jobs=n_jobs, **kwargs)
    else:
        raise NotImplementedError(f"`{metric}` has not been implemented yet.")



