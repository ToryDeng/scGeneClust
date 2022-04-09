# -*- coding: utf-8 -*-
# @Time : 2022/4/8 21:03
# @Author : Tory Deng
# @File : score.py
# @Software: PyCharm
"""
1. in-cluster scoring (how to rank genes in each cluster)
2. (optional) inter-cluster scoring (how to rank clusters)
"""
import pandas as pd


def in_cluster_score(cluster_expr: pd.DataFrame, score='var', sort: bool = True):
    if score == 'var':
        return cluster_expr.var(axis=1).sort_values(ascending=False).reset_index()
