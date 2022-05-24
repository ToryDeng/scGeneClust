# -*- coding: utf-8 -*-
# @Time : 2022/5/10 12:03
# @Author : Tory Deng
# @File : score.py
# @Software: PyCharm
import os
from multiprocessing import Pool, cpu_count
from typing import Literal

import anndata as ad
import numpy as np
import pandas as pd
from scipy.stats import kruskal
from sklearn.feature_selection import f_classif
import sys

def compute_gene_score(adata: ad.AnnData, score: Literal['f_stat', 'kw_stat']):
    if score == 'f_stat':
        adata.var[score], _ = f_classif(adata.layers['log-normalized'], adata.obs['cluster'])
        # adata.var.loc[adata.var[score] == float('inf'), score] = np.finfo(float('inf')).max-1 sys.maxsize
        adata.var.loc[adata.var[score] == float('inf'), score] = sys.maxsize
        print(float('inf') in list(adata.var[score]))
    elif score == 'kw_stat':
        cell_cluster_exprs = [adata.X[adata.obs['cluster'].values == grp, :] for grp in adata.obs['cluster'].unique()]
        with Pool(processes=cpu_count() - 1) as pool:
            result = pool.starmap(kruskal, [[expr[:, i] for expr in cell_cluster_exprs] for i in range(adata.n_vars)])
            adata.var[score] = [res[0] for res in result]  # stat, pval
    else:
        raise NotImplementedError(f"{score} has not been implemented!")
    adata.uns['gene_score'] = score
    pd.DataFrame(adata.var[score].sort_values()).to_csv('/volume1/home/yzhang/scRNA-FeatureSelection/score.csv')

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
    # pd.DataFrame(adata.var['cluster_stat'].sort_values()).to_csv('/volume1/home/yzhang/scRNA-FeatureSelection/stats.csv')
    adata.var['significant'] = adata.var['cluster_stat'] > np.quantile(cluster_stats, q=0.1)
    print('Scores of gene clusters: ', pd.Series(cluster_stats).value_counts(bins=10), sep='\n')





