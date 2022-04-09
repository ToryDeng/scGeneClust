# -*- coding: utf-8 -*-
# @Time : 2022/4/8 21:01
# @Author : Tory Deng
# @File : similarity.py
# @Software: PyCharm
"""
Metrics used to calculate the similarity between genes
"""
import pandas as pd
import numpy as np
from typing import Literal
from scipy.stats import pearsonr, spearmanr


def compute_gene_similarity(  # compute the correlation between variables
        data: pd.DataFrame,  # row: obs, col: vars
        metric: Literal['pearson', 'spearman', 'kendall']):
    if metric in ('pearson', 'spearman', 'kendall '):
        return data.corr(method=metric)
    else:
        raise NotImplementedError(f"Metric {metric} has not been implemented!")

