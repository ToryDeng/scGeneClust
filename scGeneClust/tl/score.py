# -*- coding: utf-8 -*-
# @Time : 2022/6/3 18:04
# @Author : Tory Deng
# @File : score.py
# @Software: PyCharm
import anndata as ad
import os
from multiprocessing.pool import ThreadPool
from scipy.stats import kruskal


def score_discriminative_gene(adata: ad.AnnData):
    """
    Compute the Kruskal–Wallis H test statistic or F statistic of each gene. This function is only called
    when the `mode` is 'two-way'.

    :param adata: The AnnData object
    :return: None
    """
    X, y = adata.layers['X_cell_norm'], adata.obs['cluster']
    # split expression of cell clusters
    cclusters = [X[y == grp, :] for grp in y.unique()]
    pool = ThreadPool(processes=os.cpu_count() - 1)
    result = pool.starmap(kruskal, [[expr[:, i] for expr in cclusters] for i in range(adata.n_vars)], chunksize=500)
    adata.var['score'] = [res[0] for res in result]  # stat, pval

