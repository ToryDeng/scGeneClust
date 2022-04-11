# -*- coding: utf-8 -*-
# @Time : 2022/4/8 21:01
# @Author : Tory Deng
# @File : cluster.py
# @Software: PyCharm
"""
The (hierarchical) clustering
"""

from sklearn.cluster import AgglomerativeClustering
import pandas as pd


def clustering_genes(data: pd.DataFrame, n_clusters: int):
    return AgglomerativeClustering(affinity='precomputed', n_clusters=n_clusters, linkage='average').fit_predict(data.values)
