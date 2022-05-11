# -*- coding: utf-8 -*-
# @Time : 2022/5/10 12:03
# @Author : Tory Deng
# @File : score.py
# @Software: PyCharm
import anndata as ad
import numpy as np
import pandas as pd
from typing import Literal
from sklearn.feature_selection import f_classif
from scipy.stats import kruskal
from multiprocessing import Pool


def compute_gene_score(adata: ad.AnnData, score: Literal['f_stat', 'kw_stat']):
    if score == 'f_stat':
        adata.var[score], _ = f_classif(adata.layers['log-normalized'], adata.obs['cluster'])
    elif score == 'kw_stat':
        grp_idxs = [np.argwhere(adata.obs['cluster'].values == grp) for grp in adata.obs['cluster'].unique()]
        with Pool(processes=8) as pool:
            adata.var[score] = pool.starmap(kruskal_wallis_test, [(adata.X[:, i], grp_idxs) for i in range(adata.n_vars)])
    else:
        raise NotImplementedError(f"{score} has not been implemented!")
    adata.uns['gene_score'] = score


def kruskal_wallis_test(gene_expr: np.ndarray, group_indexes: np.ndarray):
    statistic, pval = kruskal(*[gene_expr[grp_idx].squeeze() for grp_idx in group_indexes])
    return statistic


def compute_gene_cluster_score(adata: ad.AnnData, top_n: int, score: Literal['top_mean']):
    print('Evaluating gene clustering results...')
    cluster_stats = []
    for gc in adata.var['cluster'].unique():
        # find representative genes in each cluster according to gene scores
        rep_genes = adata.var.loc[adata.var['cluster'] == gc, adata.uns['gene_score']].sort_values(ascending=False).head(top_n).index
        if score == 'top_mean':
            cluster_stat = adata.var.loc[rep_genes, adata.uns['gene_score']].mean()
        else:
            raise NotImplementedError(f"{score} has not been implemented!")

        adata.var.loc[adata.var['cluster'] == gc, 'cluster_stat'] = cluster_stat
        cluster_stats.append(cluster_stat)

    adata.var['significant'] = adata.var['cluster_stat'] > np.quantile(cluster_stats, q=0.1)
    print('Scores of gene clusters: ', pd.Series(cluster_stats).value_counts(bins=10), sep='\n')





