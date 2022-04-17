# -*- coding: utf-8 -*-
# @Time : 2022/4/11 17:13
# @Author : Tory Deng
# @File : functions.py
# @Software: PyCharm
import pandas as pd

import GeneClust.steps as steps
import anndata as ad
import numpy as np
from typing import Literal, List, Union
from .utils import subset_adata


def select_genes_per_cluster(
        per_cluster_df: List[pd.DataFrame],
        n_selected_per_cluster: int
) -> np.ndarray:
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


def iteratively_select(adata: ad.AnnData, in_cluster_score, n_selected_genes, n_clusters):
    """
    Iteratively select genes from each cluster.

    :param adata: Anndata object.
    :param in_cluster_score: The type of in-cluster score of genes.
    :param n_selected_genes: The number of genes to be selected.
    :param n_clusters: The number of clusters that genes are clustered to.
    :return:
    """
    # construct the list of dataframes containing genes and their scores in each cluster, sorted by scores
    gene_info = adata.var.loc[adata.var['use_cluster'], :].copy().reset_index()
    print(f"{n_clusters - len(gene_info['cluster_label'].unique())} of all {n_clusters} clusters were dropped.")
    use_cols = ('index', in_cluster_score)
    per_cluster_gene_scores = [
        gene_info.loc[gene_info['cluster_label'] == i, use_cols].sort_values(by=in_cluster_score, ascending=False)
        for i in gene_info['cluster_label'].unique()
    ]

    # select genes from each cluster
    # first-round selection: select n_selected_genes // n_clusters genes from each cluster
    selected_genes = select_genes_per_cluster(per_cluster_gene_scores, n_selected_genes // n_clusters)
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
           dr_method: Literal['pca', 'glm-pca', 'umap'],
           n_comps: int,
           similarity: Literal['pearson', 'spearman', 'kendall', 'bayes_corr', 'mutual_info'],
           clustering: Literal['agglomerative', 'gmm'],
           n_clusters: int,
           in_cluster_score: Literal['var', 'm3drop'],
           return_genes: bool = False
           ) -> Union[pd.DataFrame, ad.AnnData]:
    """
    The main function of GeneClust.

    :param adata: Anndata object. rows represent cells, and cols represent genes.
    :param n_selected_genes: The number of genes to be selected.
    :param dr_method: The dimension reduction method.
    :param n_comps: The number of components to be used.
    :param similarity: The similarity metric.
    :param n_clusters: The number of clusters that genes are clustered to.
    :param in_cluster_score: The type of in-cluster score of genes.
    :param return_genes: Bool. whether return the selected genes (a dataframe), or the filtered anndata object.
    :return: Selected genes (an dataframe) or filtered anndata.
    """
    # dimension reduction
    adata = steps.reduce_dimension(adata, dr_method, n_comps)
    # calculate the similarity
    steps.compute_gene_similarity(adata, dr_method, similarity)
    # cluster genes
    steps.clustering_genes(adata, dr_method, similarity, clustering, n_clusters)
    # calculate in-cluster scores for genes in each cluster
    steps.in_cluster_score(adata, in_cluster_score)
    # calculate inter-cluster scores for each cluster
    steps.inter_cluster_score(adata, in_cluster_score)
    # get selected genes
    selected_genes = iteratively_select(adata, in_cluster_score, n_selected_genes, n_clusters)
    # only preserve selected genes in adata
    filtered_adata = subset_adata(adata, selected_genes, inplace=False)
    return filtered_adata.var_names.to_frame(index=False, name='Gene') if return_genes else filtered_adata
