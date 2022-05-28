# -*- coding: utf-8 -*-
# @Time : 2022/5/21 19:57
# @Author : Tory Deng
# @File : test_cluster.py
# @Software: PyCharm
import numpy as np
import pytest

from pagest.tl import reduce_dimension, cluster_genes, do_clustering
from pagest.utils import load_example_adata, set_logger


def test_cluster_genes():
    set_logger()
    adata = load_example_adata(min_cells=4)
    reduce_dimension(adata, mode='one-way', random_stat=42)
    with pytest.raises(ValueError):
        cluster_genes(adata, mode='one-way', method='gmm', n_gene_clusters=10000000, random_stat=42)
    with pytest.raises(ValueError):
        cluster_genes(adata, mode='one-way', method='leiden', n_gene_clusters=100, random_stat=42)
    cluster_genes(adata, mode='one-way', method='gmm', random_stat=42)
    cluster_genes(adata, mode='one-way', method='leiden', random_stat=42)


def test_do_clustering():
    set_logger()
    adata = load_example_adata(min_cells=2)
    reduce_dimension(adata, mode='two-way', random_stat=42)

    do_clustering(adata, 'one-way', gene_clustering_method='gmm', random_stat=42)
    assert not np.isinf(adata.var['centrality'].max())
    do_clustering(adata, 'one-way', gene_clustering_method='leiden', random_stat=42)
    assert not np.isinf(adata.var['centrality'].max())
    assert 'distances' in adata.varp and 'connectivities' in adata.varp
    do_clustering(adata, 'two-way', gene_clustering_method='gmm', random_stat=42)
    assert 'cluster' in adata.obs and 'highly_confident' in adata.obs








