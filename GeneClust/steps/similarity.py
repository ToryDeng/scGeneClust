# -*- coding: utf-8 -*-
# @Time : 2022/4/8 21:01
# @Author : Tory Deng
# @File : similarity.py
# @Software: PyCharm
"""
Metrics used to calculate the similarity between genes
"""
import pandas as pd
import anndata as ad
from typing import Literal


def compute_gene_similarity(
        adata: ad.AnnData,
        dr_method: str,
        metric: Literal['pearson', 'spearman', 'kendall']
):
    """
    Compute the similarity matrix for genes and store it in adata.varp.

    :param adata: The anndata object
    :param dr_method: The dimension reduction method
    :param metric: The similarity metric
    :return: The similarity matrix, in which each entry is the similarity between two genes
    """
    if metric in ('pearson', 'spearman', 'kendall'):
        adata.varp[metric] = adata.varm[dr_method].T.corr(method=metric).abs()  # the cols represent genes after transpose
    else:
        raise NotImplementedError(f"Metric {metric} has not been implemented!")

