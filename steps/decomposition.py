# -*- coding: utf-8 -*-
# @Time : 2022/4/8 21:00
# @Author : Tory Deng
# @File : decomposition.py
# @Software: PyCharm
"""
1. PCA
2. GLM-PCA
3. UMAP
"""
import pandas as pd
from typing import Literal
import scanpy as sc


def reduce_dimension(
        data: pd.DataFrame,  # row:  obs, col: vars
        dr_method: Literal['pca', 'glm-pca', 'umap'],
        n_componets: int
):
    if dr_method == 'pca':
        reducted_data = sc.pp.pca(data.values, n_comps=n_componets, copy=True)
        return pd.DataFrame(reducted_data, index=data.index)
    else:
        raise NotImplementedError(f"{dr_method} has not been implemented!")

