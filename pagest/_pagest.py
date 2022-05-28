# -*- coding: utf-8 -*-
# @Time : 2022/5/23 20:08
# @Author : Tory Deng
# @File : _pagest.py
# @Software: PyCharm
from typing import Literal, Optional, Union, Callable

import anndata as ad
import scanpy as sc
import numpy as np

import pagest.tl as tl
import pagest.pp as pp
from pagest.utils import set_logger


def select_from_clusters(adata: ad.AnnData, mode: Literal['one-way', 'two-way'], n_features: int):
    """
    Select `n_features` genes from gene clusters.

    :param adata: The ann oblect. adata.var must contain 'cluster' and 'centrality' (or 'stat') columns and is indexed by genes
    :param mode: `one-way` only considers patterns in genes; `two-way` considers patterns both in cells and genes
    :param n_features: The number of features to be selected
    :return: ndarray, selected features
    """
    gene_score = 'centrality' if mode == 'one-way' else 'stat'
    assert 'cluster' in adata.var, KeyError(f"Column not found: 'cluster'")
    assert gene_score in adata.var, KeyError(f"Column not found: '{gene_score}'")
    df = adata.var.loc[:, ('cluster', gene_score)]
    df.rename_axis('gene', inplace=True)

    def cluster_top_k(k: int) -> np.ndarray:
        """Select top k genes in each gene cluster."""
        return df.groupby(by='cluster')[gene_score].nlargest(k).reset_index(level=1)['gene'].values

    selected_genes = cluster_top_k(n_features // df['cluster'].unique().shape[0])
    df.drop(selected_genes, inplace=True)
    n_remained_genes = n_features - selected_genes.shape[0]
    if n_remained_genes == 0:
        return selected_genes
    else:
        while True:
            remained_genes = cluster_top_k(1)
            if n_remained_genes - remained_genes.shape[0] > 0:
                n_remained_genes -= remained_genes.shape[0]
                selected_genes = np.hstack([remained_genes, selected_genes])
                df.drop(remained_genes, inplace=True)
            elif n_remained_genes - remained_genes.shape[0] == 0:
                selected_genes = np.hstack([remained_genes, selected_genes])
                break
            else:
                randomly_chosen_genes = np.random.choice(remained_genes, size=n_remained_genes, replace=False)
                selected_genes = np.hstack([randomly_chosen_genes, selected_genes])
                break
    return selected_genes


def pagest(
        adata: ad.AnnData,
        n_features: int,
        mode: Literal['one-way', 'two-way'],
        n_components: int = 50,
        gene_clustering: Literal['gmm', 'leiden'] = 'leiden',
        n_gene_clusters: Optional[int] = None,
        gene_distance: Union[str, Callable] = None,
        n_gene_neighbors: Optional[int] = 30,
        n_cell_clusters: Optional[int] = None,
        gene_cluster_score: Literal['silhouette', 'spearman'] = 'spearman',
        drop_quantile: float = 0.1,
        stat: Literal['kw', 'f'] = 'kw',
        subset: bool = False,
        verbose: Literal[0, 1, 2] = 1,
        random_stat: Optional[int] = None
) -> Optional[np.ndarray]:
    """
    Pattern-aware Gene Selection (PAGEST). If `mode='one-way'`, it will only find patterns in genes and preserve the
    most representative genes to reduce the redundancy in genes. If `mode='two-way'`, it will simultaneously find patterns
    in cells and genes, and select genes that are not only pattern-specific but also discriminative.

    :param adata: The AnnData object
    :param n_features: The number of features to be selected
    :param mode: `one-way` only considers patterns in genes; `two-way` considers patterns both in cells and genes
    :param n_components: The number of used principle components
    :param gene_clustering: The gene clustering method
    :param n_gene_clusters: The number of gene clusters. Will be automatically determined if None.
    :param gene_distance: The metric used to measure the distances between genes
    :param n_gene_neighbors: The number of gene neighbors. Only used in leiden clustering.
    :param n_cell_clusters: The number of cell clusters
    :param gene_cluster_score: 'silhouette': mean silhouette scores; 'spearman': mean absolute spearman correlation
    :param drop_quantile: The quantile of number of gene clusters to compute
    :param stat: The kind of statistic used to represent differences of gene expression between cell groups
    :param subset: Inplace subset to informative genes if True otherwise merely return the selected genes
    :param verbose: An integer in [0, 1, 2]. 0: only print warnings and errors; 1: also print info; 2: also print debug message
    :param random_stat: The random seed. Default is none
    :return: (Optional) the selected genes if `subset` is False
    """
    np.random.seed(random_stat)
    set_logger(verbose)

    pp.preprocess(adata)
    pp.reduce_dimension(adata, mode, n_components, random_stat)
    tl.do_clustering(adata, mode, gene_clustering, n_gene_clusters, gene_distance, n_gene_neighbors, n_cell_clusters, random_stat)
    # tl.score_gene_cluster(adata, gene_cluster_score)
    filtered_adata = tl.filter_adata2(adata, mode, drop_quantile)

    if mode == 'two-way':
        tl.score_discriminative_gene(filtered_adata, stat)
        adata.var['stat'] = filtered_adata.var['stat']
    selected_genes = select_from_clusters(filtered_adata, mode, n_features)
    is_informative = np.isin(adata.var_names, selected_genes)
    if is_informative.sum() != n_features:
        raise RuntimeError(f"Only found {is_informative.sum()} informative genes in adata.var_names, not {n_features}. "
                           f"Not used gene(s): {np.isin(selected_genes, adata.var_names).sum()}")

    if subset:
        adata._inplace_subset_var(is_informative)
    else:
        return selected_genes











