# -*- coding: utf-8 -*-
# @Time : 2022/4/11 17:13
# @Author : Tory Deng
# @File : functions.py
# @Software: PyCharm
from typing import Literal, List, Union

import anndata as ad
import numpy as np
import pandas as pd

import GeneClust.steps as steps
from .utils import subset_adata


def select_genes_per_cluster(per_cluster_df: List[pd.DataFrame], n_selected_per_cluster: int) -> np.ndarray:
    """
    Select top genes from each cluster.

    :param per_cluster_df: A list of sorted dataframe. In each dataframe, the first column contains gene names, and
    the second columns contains gene scores in the cluster.
    :param n_selected_per_cluster: The number of genes to be selected from each cluster.
    :return: 1-D array. All genes selected from clusters.
    """
    per_cluster_selected_genes = []
    for i, cluster_gene_scores in enumerate(per_cluster_df.copy()):
        if cluster_gene_scores.shape[0] >= n_selected_per_cluster:
            top_selected_genes = cluster_gene_scores.iloc[:n_selected_per_cluster, 0]
            per_cluster_selected_genes.append(top_selected_genes.values)
            per_cluster_df[i].drop(top_selected_genes.index, inplace=True)
        else:
            per_cluster_selected_genes.append(cluster_gene_scores['index'].values)
            per_cluster_df[i].drop(cluster_gene_scores['index'].index, inplace=True)
    return np.concatenate(per_cluster_selected_genes)


def iteratively_select(adata: ad.AnnData, n_selected_genes) -> np.ndarray:
    """
    Iteratively select genes from each cluster.

    :param adata: Anndata object.
    :param n_selected_genes: The number of genes to be selected.
    :return:
    """
    # construct the list of dataframes containing genes and their scores in each cluster, sorted by scores
    n_total_clusters = adata.var['cluster_label'].unique().shape[0]
    preserved_genes = adata.var.loc[adata.var['use_cluster'], :].copy().reset_index()
    preserved_clusters = preserved_genes['cluster_label'].unique()

    print(f"{n_total_clusters - preserved_clusters.shape[0]} of all {n_total_clusters} clusters were dropped.")
    use_cols = ('index', adata.uns['intra_cluster_score'])
    per_cluster_gene_scores = [
        preserved_genes.loc[preserved_genes['cluster_label'] == i, use_cols].sort_values(
            by=adata.uns['intra_cluster_score'], ascending=False
        ) for i in preserved_clusters
    ]

    # select genes from each cluster
    # first-round selection: select n_selected_genes // n_clusters genes from each cluster
    selected_genes = select_genes_per_cluster(per_cluster_gene_scores, n_selected_genes // n_total_clusters)
    # if lack some genes, keep selecting from each cluster
    while selected_genes.shape[0] < n_selected_genes:
        single_selection = select_genes_per_cluster(per_cluster_gene_scores, 1)
        if single_selection.shape[0] > n_selected_genes - selected_genes.shape[0]:
            choice_idxs = np.random.choice(single_selection.shape[0], n_selected_genes - selected_genes.shape[0], False)
            single_selection = single_selection[choice_idxs]
        selected_genes = np.append(selected_genes, single_selection)

    return selected_genes


def select(adata: ad.AnnData,
           n_selected_genes: int,
           dr_method: Literal['pca', 'umap', 'pca-umap'],
           n_comps: int,
           distance: Literal['pearson', 'spearman', 'kendall', 'bayesian', 'euclidean', 'mahalanobis', 'rho_p', 'phi_s'],
           clustering: Literal['agg', 'gmm'],
           n_clusters: int,
           in_cluster_score: Literal['m3drop', 'seurat', 'center'],
           inter_cluster_score: Literal['top3', 'silhouette'],
           return_genes: bool = False
           ) -> Union[pd.DataFrame, ad.AnnData]:
    """
    The main function of GeneClust.

    :param adata: Anndata object. rows represent cells, and cols represent genes.
    :param n_selected_genes: The number of genes to be selected.
    :param dr_method: The dimension reduction method.
    :param n_comps: The number of components to be used.
    :param distance: The distance metric.
    :param clustering: The clustering method.
    :param n_clusters: The number of clusters that genes are clustered to.
    :param in_cluster_score: The type of in-cluster score of genes in each cluster.
    :param inter_cluster_score: The type of inter-cluster score of each cluster.
    :param return_genes: Bool. Whether return the selected genes (a dataframe), or the filtered anndata object.
    :return: Selected genes (an dataframe) or filtered anndata.
    """
    # dimension reduction
    adata = steps.reduce_dimension(adata, dr_method, n_comps)
    # calculate the distance. when using gmm and the inter_cluster_score is not silhouette score, we don't compute
    if clustering != 'gmm' or inter_cluster_score == 'silhouette':
        steps.compute_gene_distance(adata, distance)
    # cluster genes
    steps.clustering_genes(adata, clustering, n_clusters)
    # calculate in-cluster scores for genes in each cluster
    steps.in_cluster_score(adata, in_cluster_score)
    # calculate inter-cluster scores for each cluster
    steps.inter_cluster_score(adata, inter_cluster_score)
    # get selected genes
    selected_genes = iteratively_select(adata, n_selected_genes)
    # only preserve selected genes in adata
    filtered_adata = subset_adata(adata, selected_genes, inplace=False)
    return filtered_adata.var_names.to_frame(index=False, name='Gene') if return_genes else filtered_adata
