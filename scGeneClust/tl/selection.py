# -*- coding: utf-8 -*-
# @Time : 2023/1/6 23:52
# @Author : Tory Deng
# @File : selection.py
# @Software: PyCharm
from typing import Literal

import anndata as ad
import numpy as np
from loguru import logger
from sklearn.ensemble import IsolationForest


def select_from_clusters(
        adata: ad.AnnData,
        version: Literal['fast', 'ps'],
        post_hoc_filtering: bool = True,
        random_state: int = 0
):
    """
    If `version='fast'`, select the farthest and closest genes in each non-singleton clusters.
    If `version='ps'`, select the representative genes in each high-density clusters.

    Parameters
    ----------
    adata : AnnData
        The annotated data matrix of shape `n_obs` × `n_vars`.
        Rows correspond to cells and columns to genes.
    version : str
        Choose the version of GeneClust.
    post_hoc_filtering : bool
        Whether to find outliers in singleton gene clusters (GeneClust-fast) or low-density genes (GeneClust-ps)
        after gene clustering.
    random_state : int
        Change to use different initial states for the optimization.

    Returns
    -------
    selected_features: ndarray
        Names of selected genes.
    """
    gene_cluster_counts = adata.var['cluster'].value_counts()
    logger.debug(
        f"Gene cluster size: \n{gene_cluster_counts.rename('size').to_frame().rename_axis('cluster', axis=0)}"
    )
    if version == 'fast':
        is_single_cluster = adata.var['cluster'].isin(gene_cluster_counts[gene_cluster_counts == 1].index)  # (n_genes,)
        if post_hoc_filtering and (n_single_clusters := is_single_cluster.sum()) > 0:  # if contains single clusters
            single_cluster_genes = adata.var_names[is_single_cluster]
            deviances = compute_deviance(adata.X[:, is_single_cluster])  # (n_cells, n_single_genes)
            is_outlier = IsolationForest(random_state=random_state).fit_predict(deviances.reshape(-1, 1)) == -1  # (n_single_genes,)
            logger.opt(colors=True).debug(f"Found <yellow>{is_outlier.sum()}</yellow> outlier(s) in <yellow>{n_single_clusters}</yellow> single gene cluster(s).")
            outlier_genes = single_cluster_genes[is_outlier]
        else:
            outlier_genes = np.array([])
        # select genes from each cluster
        grouped = adata.var.loc[~is_single_cluster, ('cluster', 'closeness')].rename_axis('gene').groupby(by='cluster')['closeness']
        max_genes = grouped.nlargest(1).reset_index(level=1)['gene'].values
        min_genes = grouped.nsmallest(1).reset_index(level=1)['gene'].values
        selected_features = np.hstack((max_genes, min_genes, outlier_genes))
    else:
        is_low_density_gene = adata.var['cluster'] == -1  # (n_genes,)
        if (n_low_density_clusters := is_low_density_gene.sum()) > 0:
            logger.opt(colors=True).debug(f"Found <yellow>{n_low_density_clusters}</yellow> low density genes.")
            adata.var['representative'] = False
            grouped = adata.var[~is_low_density_gene].groupby(by='cluster')['relevance']
            largest_genes = grouped.nlargest(1).index.get_level_values(level=1)
            adata.var.loc[largest_genes, 'representative'] = True
            logger.opt(colors=True).debug(
                f"<yellow>{adata.var['representative'].sum()}</yellow> representative gene(s) in high-density gene cluster(s)."
            )
            if post_hoc_filtering:
                high_density_values = adata.var.loc[~is_low_density_gene, 'outlier_score']
                threshold = high_density_values[high_density_values > 0].median()
                is_outlier = adata.var.loc[is_low_density_gene, 'outlier_score'] >= threshold
                adata.var.loc[is_low_density_gene, 'representative'] = is_outlier
                logger.opt(colors=True).debug(f"Found <yellow>{is_outlier.sum()}</yellow> outliers in low-density gene cluster(s).")
        selected_features = adata.var_names[adata.var['representative']].to_numpy()
    return selected_features


def compute_deviance(X: np.ndarray):
    """Compute deviance score for each single gene cluster."""
    pi = X.sum(0) / X.sum()
    n = X.sum(1)[:, np.newaxis]
    with np.errstate(all='ignore'):
        left_half = np.nansum(X * np.log(X / n.dot(pi[np.newaxis, :])), axis=0)
        right_half = np.nansum((n - X) * (np.log((n - X) / (n.dot((1 - pi)[np.newaxis, :])))), axis=0)
    return 2 * (left_half + right_half)

