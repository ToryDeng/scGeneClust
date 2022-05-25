# -*- coding: utf-8 -*-
# @Time : 2022/5/23 11:52
# @Author : Tory Deng
# @File : test_score.py
# @Software: PyCharm
from pagest.tl import reduce_dimension, do_clustering, score_gene_cluster
from pagest.utils import load_example_adata, set_logger


def test_score_gene_cluster():
    set_logger()
    adata = load_example_adata(min_cells=4)
    reduce_dimension(adata, mode='two-way', random_stat=42)
    do_clustering(adata, 'two-way', gene_clustering_method='gmm', random_stat=42)

    score_gene_cluster(adata, method='silhouette')
    score_gene_cluster(adata, method='spearman')
