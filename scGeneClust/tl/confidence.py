# -*- coding: utf-8 -*-
# @Time : 2023/1/5 22:19
# @Author : Tory Deng
# @File : confidence.py
# @Software: PyCharm
import contextlib
import os
import random
import torch
from functools import partial
from multiprocessing.pool import ThreadPool
# from SpaGCN import SpaGCN, calculate_adj_matrix, search_l, search_res
from typing import Literal
from typing import Optional

import SpaGCN as spg
import anndata as ad
import igraph as ig
import leidenalg
import numpy as np
import pandas as pd
import squidpy as sq
from loguru import logger
from scipy.stats import entropy
from sklearn.mixture import GaussianMixture

X_pca = None


def find_high_confidence_cells(
        adata: ad.AnnData,
        n_cell_clusters: int,
        n_components: int = 10,
        max_workers: int = os.cpu_count() - 1,
        random_state: int = 0
):
    """
    Find cells which belong to certain cluster with high confidence.

    Parameters
    ----------
    adata
        The annotated data matrix of shape `n_obs` × `n_vars`.
        Rows correspond to cells and columns to genes.
    n_cell_clusters
        The number of clusters in cell clustering used to find high-confidence cells. Only valid in GeneClust-ps.
    n_components
        The number of principal components used along with the first component. Only valid in GeneClust-ps.
    max_workers
        The maximum value of workers which can be used during feature selection.
    random_state
        Change to use different initial states for the optimization.
    """
    logger.info(f"Finding high-confidence cells...")
    # compute the frequency matrix
    global X_pca
    X_pca = adata.obsm['X_pca']
    partial_compute = partial(_compute_cell_co_membership, n_clusters=n_cell_clusters, random_state=random_state)
    with ThreadPool(processes=max_workers) as pool:
        results = pool.map(partial_compute, range(2, 2 + n_components))
    frequency_matrix = np.sum(results, axis=0)  # nonzero values only exist in the upper right part
    # find threshold
    for freq_th in range(8, 0, -1):
        cut_matrix = np.where(frequency_matrix < freq_th, 0, frequency_matrix)
        cluster_labels = leiden(cut_matrix, seed=random_state)
        cluster_counts = pd.Series(cluster_labels).value_counts(ascending=False)
        if (cluster_counts < 10).any():  # has_small_cluster
            cut_k = min(n_cell_clusters, np.argwhere((cluster_counts < 10).values).squeeze(axis=1)[0])
        else:
            cut_k = cluster_counts.shape[0]
        is_confident = np.isin(cluster_labels, cluster_counts.index[:cut_k])
        if is_confident.sum() / adata.n_obs > 0.1:
            logger.opt(colors=True).debug(f"Final frequency cutoff: <yellow>{freq_th}</yellow>")
            break
    adata._inplace_subset_obs(is_confident)
    adata.obs['cluster'] = cluster_labels.astype(str)[is_confident]
    logger.opt(colors=True).info(
        f"Found <yellow>{adata.n_obs}</yellow> (<yellow>{np.round(adata.n_obs / is_confident.shape[0] * 100)}%</yellow>) high-confidence cells."
    )


def _compute_cell_co_membership(idx: np.ndarray, n_clusters: int, random_state: int, p: float = 0.95) -> np.ndarray:
    """
    Perform GMM clustering on certain PCs.

    Parameters
    ----------
    idx : np.ndarray
        Number of gene-level principal components used in GMM
    n_clusters : int
        Number of cell clusters
    random_state : int
        Change to use different initial states for the optimization.
    p : float, default=0.95
        The probability threshold for a certain cell to be considered to belong to the cluster.

    Returns
    -------
    co_membership : ndarray
        Contains only 0 (not belong to same cluster) and 1 (belong to same cluster) with shape (n_cells, n_cells).
    """
    # prepare data
    global X_pca
    X = X_pca[:, :idx]
    if X.ndim == 1:
        X = X.reshape(-1, 1)
    # GMM clustering
    gmm = GaussianMixture(n_components=n_clusters, init_params='k-means++', random_state=random_state)
    gmm.fit(X)
    cell_labels, cell_probas = gmm.predict(X), gmm.predict_proba(X).max(1)
    co_membership = np.zeros((X.shape[0], X.shape[0]))  # rows correspond to cells
    for i in range(X.shape[0] - 1):
        if cell_probas[i] >= p:
            co_membership[i, i + 1:] = np.logical_and(cell_probas[i + 1:] > p, cell_labels[i + 1:] == cell_labels[i])
    return co_membership


def leiden(adjacency: np.ndarray, resolution: float = 1.0, seed: Optional[int] = None):
    """
    Create an undirected graph and perform leiden clustering on a given adjacency matrix.

    Parameters
    ----------
    adjacency : np.ndarray
        An adjacency matrix. Nonzero values only exist in the upper right part.
    resolution : float
        The resolution parameter.
    seed : Optional[int]
        Seed for the random number generator.

    Returns
    -------
    cluster_labels
        Cluster label of each cell.
    """
    G = ig.Graph.Weighted_Adjacency(adjacency, mode="upper")  # igraph checks adjacency matrix internally
    logger.debug("Leiden clustering starts...")
    partition = leidenalg.find_partition(
        G,
        partition_type=leidenalg.RBConfigurationVertexPartition,
        weights=G.es['weight'],
        n_iterations=-1,
        resolution_parameter=resolution,
        seed=seed
    )
    cluster_labels = np.array(partition.membership)
    logger.debug("Leiden clustering finished!")
    return cluster_labels


def run_spaGCN(
        adata: ad.AnnData,
        img: np.ndarray,
        n_spot_cluster: int,
        shape: Literal['hexagon', 'square'] = 'hexagon',
        random_state: int = 0
) -> np.ndarray:
    logger.opt(colors=True).debug(f"spaGCN starts running on <yellow>{adata.n_obs}</yellow> spots and <yellow>{adata.n_vars}</yellow> genes...")
    # prepare positional information
    x_array, y_array = adata.obs["array_row"].values, adata.obs["array_col"].values
    x_pixel, y_pixel = adata.obsm['spatial'][:, 1], adata.obsm['spatial'][:, 0]

    with open(os.devnull, "w") as f, contextlib.redirect_stdout(f):
        if img is None:
            adj = spg.calculate_adj_matrix(x=x_pixel, y=y_pixel, histology=False)
        else:
            adj = spg.calculate_adj_matrix(
                x=x_pixel, y=y_pixel, x_pixel=x_pixel, y_pixel=y_pixel, image=img, beta=49, alpha=1, histology=True
            )
        l = spg.search_l(0.5, adj, start=0.01, end=1000, tol=0.01, max_run=100)
        res = spg.search_res(adata, adj, l, n_spot_cluster, start=0.7, step=0.1, tol=5e-3, lr=0.05, max_epochs=20,
                         r_seed=random_state, t_seed=random_state, n_seed=random_state)
        clf = spg.SpaGCN()
        clf.set_l(l)
        random.seed(random_state)
        torch.manual_seed(random_state)
        np.random.seed(random_state)
        clf.train(adata, adj, init_spa=True, init="louvain", res=res, tol=5e-3, lr=0.05, max_epochs=200)
        y_pred, prob = clf.predict()
        cluster_labels = spg.spatial_domains_refinement_ez_mode(
            sample_id=adata.obs.index, pred=y_pred, x_array=x_array, y_array=y_array, shape=shape
        )
    logger.opt(colors=True).debug("<magenta>spaGCN</magenta> completed.")
    return np.array(cluster_labels)


def find_high_confidence_spots(
        adata: ad.AnnData,
        img: np.ndarray,
        n_spot_cluster: int,
        shape: Literal['hexagon', 'square'] = 'hexagon',
        n_rings: int = 2,
        random_state: int = 0
):
    logger.info(f"Finding high-confidence spots...")
    adata.obs["cluster"] = run_spaGCN(adata, img, n_spot_cluster, shape, random_state)

    if shape == 'hexagon':
        n_neighs, min_neighs = 6, 1.5 * n_rings * (n_rings + 1)
    else:
        n_neighs, min_neighs = 4, 1.0 * n_rings * (n_rings + 1)
    sq.gr.spatial_neighbors(adata, n_rings=n_rings, coord_type="grid", n_neighs=n_neighs)
    n_true_neighbors = adata.obsp["spatial_connectivities"].sum(0).A1
    neigh_entropies, is_same_cluster = [], []
    spots_clusters = adata.obs['cluster'].values
    for i in range(adata.n_obs):
        neigh_idx = adata.obsp["spatial_connectivities"].getrow(i).nonzero()[1]
        if neigh_idx.shape[0] == 0:  # th spot doesn't have any neighbor
            is_same_cluster.append(False)
            neigh_entropies.append(np.inf)
            continue
        unique_clusters, counts = np.unique(spots_clusters[neigh_idx], return_counts=True)  # count the neighborhood
        neigh_main_clusters = unique_clusters[np.argmax(counts)].flatten()
        if spots_clusters[i] in neigh_main_clusters:
            is_same_cluster.append(True)
        else:
            is_same_cluster.append(False)
        neigh_entropies.append(entropy(counts))

    neigh_entropies, is_same_cluster = np.array(neigh_entropies), np.array(is_same_cluster)
    entropy_threshold = np.percentile(neigh_entropies, 20)
    logger.opt(colors=True).debug(f"Entropy threshold: <yellow>{entropy_threshold}</yellow>")
    is_confident = (n_true_neighbors >= min_neighs) & is_same_cluster & (neigh_entropies <= entropy_threshold)
    hc_clusters, hc_cluster_counts = np.unique(spots_clusters[is_confident], return_counts=True)
    small_hc_clusters = hc_clusters[hc_cluster_counts < 10]

    adata._inplace_subset_obs(is_confident & (~adata.obs['cluster'].isin(small_hc_clusters)))
    logger.opt(colors=True).info(
        f"Found <yellow>{adata.n_obs}</yellow> (<yellow>{np.round(adata.n_obs / is_confident.shape[0] * 100)}%</yellow>) high-confidence spots."
    )

