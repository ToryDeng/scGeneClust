# -*- coding: utf-8 -*-
# @Time : 2022/6/3 17:09
# @Author : Tory Deng
# @File : _utils.py
# @Software: PyCharm
import anndata as ad
import numpy as np
import leidenalg
import scanpy as sc
from loguru import logger
from typing import Literal, Union, Optional
from scGeneClust.tl.distance import compute_pairwise_distances
from sklearn.metrics.pairwise import paired_distances
from sklearn.preprocessing import minmax_scale
from sklearn.neighbors import kneighbors_graph
from sklearn.mixture import GaussianMixture


def _compute_1w_gene_score(adata: ad.AnnData, centers: np.ndarray) -> np.ndarray:
    """
    Compute the distance of each gene to its cluster mean, min-max normalized in cluster.
    The gene scores is 1 - normalized distances.

    :param adata: The AnnData object
    :param centers: The cluster means
    :return: distances of all genes
    """
    all_distances = paired_distances(adata.varm['pca'], centers[adata.var['cluster']])
    for gene_cluster in np.unique(adata.var['cluster']):
        gene_cluster_mask = adata.var['cluster'] == gene_cluster
        all_distances[gene_cluster_mask] = 1 - minmax_scale(all_distances[gene_cluster_mask])
    return all_distances


def _find_neighbors(
        adata: ad.AnnData,
        on: Literal['cell', 'gene'],
        dis_metric: str,
        n_neighbors: int,
):
    attr = 'obs' if on == 'cell' else 'var'
    attrp, attrm = attr + 'p', attr + 'm'

    logger.debug(f"Computing {on} distances using {dis_metric}...")

    getattr(adata, attrp)['distances'] = compute_pairwise_distances(getattr(adata, attrm)['pca'], metric=dis_metric)
    logger.debug(f"Finding KNNs on {on}...")
    knn_graph = kneighbors_graph(
        getattr(adata, attrp)['distances'], n_neighbors, metric='precomputed', n_jobs=-1
    ).toarray()
    logger.debug(f"Finding MNNs on {on}...")
    mnn_graph = np.logical_and(knn_graph, knn_graph.T).astype(int)
    getattr(adata, attrp)['connectivities'] = mnn_graph
    logger.debug(f"Neighbor finding finished!")


def leiden(
        data: Union[ad.AnnData, np.ndarray],
        on: Literal['cell', 'gene'] = None,
        resolution: float = 1.0,
        seed: Optional[int] = None
):
    if isinstance(data, ad.AnnData):
        if on is None:
            raise TypeError("`on` can not be None when `data` is `ad.AnnData` object!")
        attr = 'obs' if on == 'cell' else 'var'
        attrp, attrm = attr + 'p', attr + 'm'
        if 'connectivities' not in getattr(data, attrp):
            raise ValueError(f"Can not find `connectivities` in {attrp}. Please calculate {on} neighbors first.")
        adjacency = getattr(data, attrp)['connectivities']
    elif isinstance(data, np.ndarray):
        if data.ndim != 2:
            raise ValueError(f"Dimension of data must be 2, but got {data.ndim} dim(s) array!")
        if data.shape[0] != data.shape[1]:
            raise ValueError(f"`data` should be a square matrix!")
        adjacency = data
    else:
        raise TypeError(f"`data` can only be `ad.AnnData` or `np.ndarray`!")

    directed = False if np.allclose(adjacency, adjacency.T) else True
    logger.debug(f"Creating {'directed' if directed else 'undirected'} graph on {on}s...")
    G = sc._utils.get_igraph_from_adjacency(adjacency, directed=directed)
    logger.debug("Leiden clustering starts...")
    partition = leidenalg.find_partition(G, partition_type=leidenalg.RBConfigurationVertexPartition,
                                         weights=G.es['weight'],
                                         n_iterations=-1, resolution_parameter=resolution, seed=seed)
    cluster_labels = np.array(partition.membership)
    logger.debug("Leiden clustering finished!")
    return cluster_labels


def _compute_cell_co_membership(
        X: np.ndarray,
        n_clusters: int,
        seed: Optional[int] = None,
        p: float = 0.95
) -> np.ndarray:
    X = X.reshape(-1, 1) if X.ndim == 1 else X
    gmm = GaussianMixture(n_components=n_clusters, init_params='k-means++', random_state=seed)
    gmm.fit(X)
    cluster_labels = gmm.predict(X)
    dense_co_membership = []
    cell_proba = gmm.predict_proba(X)
    max_proba = cell_proba.max(1)
    for cell_i in range(X.shape[0] - 1):
        for cell_j in range(cell_i + 1, X.shape[0]):
            if cluster_labels[cell_i] == cluster_labels[cell_j] and max_proba[cell_i] > p and max_proba[cell_j] > p:
                dense_co_membership.append(1.)
            else:
                dense_co_membership.append(0.)
    return np.array(dense_co_membership)
