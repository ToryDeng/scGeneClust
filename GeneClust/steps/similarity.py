# -*- coding: utf-8 -*-
# @Time : 2022/4/8 21:01
# @Author : Tory Deng
# @File : similarity.py
# @Software: PyCharm
"""
Metrics used to calculate the similarity between genes
"""
import pandas as pd
from typing import Literal


def compute_gene_similarity(
        data: pd.DataFrame,
        metric: Literal['pearson', 'spearman', 'kendall']
) -> pd.DataFrame:
    """
    Compute the similarity matrix for genes.

    :param data: The latent representations of genes (columns)
    :param metric: The similarity metric
    :return: The similarity matrix, in which each entry is the similarity between two genes
    """
    if metric in ('pearson', 'spearman', 'kendall'):
        return data.corr(method=metric)
    else:
        raise NotImplementedError(f"Metric {metric} has not been implemented!")

