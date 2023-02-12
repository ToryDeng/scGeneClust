# -*- coding: utf-8 -*-
# @Time : 2023/1/6 2:19
# @Author : Tory Deng
# @File : information.py
# @Software: PyCharm
import os
from itertools import combinations
from multiprocessing.pool import Pool
from typing import Tuple

import anndata as ad
import igraph as ig
import numpy as np
from loguru import logger
from scipy.spatial.distance import squareform
from sklearn.feature_selection import mutual_info_classif, mutual_info_regression

expr_mtx, clusters, seed = None, None, 0


def _compute_relevance(i: int):
    global expr_mtx, clusters, seed
    g1 = expr_mtx[:, i].reshape(-1, 1)
    return mutual_info_classif(g1, clusters, discrete_features=False, random_state=seed)[0]


def find_relevant_genes(adata: ad.AnnData, top_pct: int, max_workers: int = os.cpu_count() - 1, random_state: int = 0):
    """
    Compute relevance of each gene to the pseudo cell types (stored in `adata.var['relevance']`),
    and then inplace subsets genes with the highest relevance as relevant genes.

    Parameters
    ----------
    adata : AnnData
        The annotated data matrix of shape `n_obs` × `n_vars`.
        Rows correspond to cells and columns to genes.
    top_pct : int
        The percentage of relevant genes. This parameter should be between 0 and 100.
    max_workers : int
        The maximum value of workers which can be used during feature selection.
    random_state : int
        Change to use different initial states for the optimization.
    """
    logger.info("Finding relevant genes...")
    global expr_mtx, clusters, seed
    expr_mtx, clusters, seed = adata.X, adata.obs['cluster'], random_state
    with Pool(processes=max_workers) as pool:
        relevance = np.array(pool.map(_compute_relevance, range(adata.n_vars)))
    logger.debug("Gene relevance computed!")
    rlv_th = np.percentile(relevance, 100 - top_pct)
    logger.opt(colors=True).debug(f"Relevance threshold: <yellow>{rlv_th}</yellow>")
    is_relevant = relevance > rlv_th
    adata._inplace_subset_var(is_relevant)
    adata.var['relevance'] = relevance[is_relevant]
    logger.opt(colors=True).info(
        f"<yellow>{is_relevant.sum()}</yellow> (<yellow>{top_pct}%</yellow>) genes are marked as relevant genes."
    )


def _compute_redundancy(gene_index_pair: Tuple[int, int]):
    global expr_mtx, clusters, seed
    g1, g2 = expr_mtx[gene_index_pair[0], :].reshape(-1, 1), expr_mtx[gene_index_pair[1], :]
    return mutual_info_regression(g1, g2, discrete_features=False, random_state=seed)[0]


def compute_gene_redundancy(adata: ad.AnnData, max_workers: int = os.cpu_count() - 1, random_state: int = 0):
    """
    Compute redundancy of all relevant gene pairs and store them as a symmetric matrix in `adata.varp['redundancy']`.

    Parameters
    ----------
    adata : AnnData
        The annotated data matrix of shape `n_obs` × `n_vars`.
        Rows correspond to cells and columns to genes.
    max_workers : int
        The maximum value of workers which can be used during feature selection.
    random_state : int
        Change to use different initial states for the optimization.
    """
    logger.info("Computing gene redundancy...")
    global expr_mtx, clusters, seed
    expr_mtx, clusters, seed = adata.varm['X_pca'], None, random_state
    with Pool(processes=max_workers) as pool:  # upper triangular matrix
        adata.varp['redundancy'] = squareform(pool.map(_compute_redundancy, combinations(range(adata.n_vars), 2)))
    logger.info(f"Gene redundancy computed.")


def _compute_complementarity(gene_index_pair: Tuple[int, int]):
    global expr_mtx, clusters, seed
    i, j = gene_index_pair[0], gene_index_pair[1]

    cmi = 0
    for clus in np.unique(clusters):
        clus_mask = clusters == clus
        g1, g2 = expr_mtx[clus_mask, i].reshape(-1, 1), expr_mtx[clus_mask, j]
        rlv = mutual_info_regression(g1, g2, discrete_features=False, n_neighbors=3, random_state=seed)[0]
        cmi += rlv * clus_mask.sum() / expr_mtx.shape[0]
    return cmi


def compute_gene_complementarity(adata: ad.AnnData, max_workers: int = os.cpu_count() - 1, random_state: int = 0):
    """
    Build the gene-gene redundancy graph first, and generate an MST from the graph.
    Then compute complementarity of any two genes that are connected in the MST.
    Store the pairwise gene complementarity in `adata.uns['mst_edges_complm']` as an ndarray of shape (n_edges, ).
    Store the edges of MST in `adata.uns['mst_edges']` as an ndarray of shape (n_edges, 2).

    Parameters
    ----------
    adata : AnnData
        The annotated data matrix of shape `n_obs` × `n_vars`.
        Rows correspond to cells and columns to genes.
    max_workers : int
        The maximum value of workers which can be used during feature selection.
    random_state : int
        Change to use different initial states for the optimization.
    """
    logger.info("Computing gene complementarity...")
    global expr_mtx, clusters, seed
    expr_mtx, clusters, seed = adata.X, adata.obs['cluster'], random_state
    #  build MST first
    logger.debug("Building MST...")

    G = ig.Graph.Weighted_Adjacency(-adata.varp['redundancy'], mode="undirected")
    MST = G.spanning_tree(weights=G.es["weight"])
    adata.uns['mst_edges'] = MST.get_edge_dataframe()[['source', 'target']].values

    logger.debug("MST built.")
    # compute complementarity in MST
    with Pool(processes=max_workers) as pool:
        adata.uns['mst_edges_complm'] = np.array(pool.map(_compute_complementarity, MST.get_edgelist()))
    logger.info(f"Gene complementarity computed.")
