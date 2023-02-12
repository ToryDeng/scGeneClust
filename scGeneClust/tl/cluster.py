# -*- coding: utf-8 -*-
# @Time : 2023/1/5 8:34
# @Author : Tory Deng
# @File : cluster.py
# @Software: PyCharm
import os
from typing import Literal

import anndata as ad
import numpy as np
from hdbscan._hdbscan_linkage import label
from hdbscan._hdbscan_tree import condense_tree, compute_stability, get_clusters, outlier_scores
from loguru import logger
from sklearn.cluster import MiniBatchKMeans
from sklearn.metrics.pairwise import paired_distances
from sklearn.preprocessing import minmax_scale

from .confidence import find_high_confidence_cells, find_high_confidence_spots
from .information import find_relevant_genes, compute_gene_redundancy, compute_gene_complementarity


def cluster_genes(
        adata: ad.AnnData,
        img: np.ndarray,
        version: Literal['fast', 'ps'],
        modality: Literal['sc', 'st'] = 'sc',
        shape: Literal['hexagon', 'square'] = 'hexagon',
        n_gene_clusters: int = None,
        n_obs_clusters: int = None,
        n_components: int = 10,
        relevant_gene_pct: int = 20,
        max_workers: int = os.cpu_count() - 1,
        random_state: int = 0
):
    """
    Cluster genes using mini-batch k-means (GeneClust-fast) or an graph-based algorithm (GeneClust-ps).

    Parameters
    ----------
    adata : AnnData
        The annotated data matrix of shape `n_obs` × `n_vars`.
        Rows correspond to cells and columns to genes.
    img : ndarray
        The image of tissue section.
    version : str
        Choose the version of GeneClust.
    modality : Literal['sc', 'st'], default='sc'
        Type of the dataset. 'sc' for scRNA-seq data, 'st' for spatially resolved transcriptomics (SRT) data.
    shape : Literal['hexagon', 'square'], default='hexagon'
        The shape of spot neighbors. 'hexagon' for Visium data, 'square' for ST data.
    n_gene_clusters : int
        The number of gene clusters used to cluster genes. Only valid in GeneClust-fast.
    n_obs_clusters : int
        The number of cell clusters used to find high-confidence cells. Only valid in GeneClust-ps.
    n_components : int, default=10
        The number of principal components used along with the first component. Only valid in GeneClust-ps.
    relevant_gene_pct: int, default=20
        The percentage of relevant genes. This parameter should be between 0 and 100. Only valid in GeneClust-ps.
    max_workers: int, default=os.cpu_count() - 1
        The maximum value of workers which can be used during feature selection.
    random_state : int, default=0
        Change to use different initial states for the optimization.
    """
    logger.info("Clustering genes...")
    if version == 'fast':
        km = MiniBatchKMeans(
            n_clusters=n_gene_clusters, batch_size=max(1024, adata.n_vars // 10), random_state=random_state,
            n_init='auto'
        )
        adata.var['cluster'] = km.fit_predict(adata.varm['X_pca'])
        adata.var['closeness'] = compute_gene_closeness(adata, km.cluster_centers_)
    else:
        if modality == 'sc':
            find_high_confidence_cells(adata, n_obs_clusters, n_components, max_workers, random_state)
        else:
            find_high_confidence_spots(adata, img, n_obs_clusters, shape, random_state=random_state)
        find_relevant_genes(adata, relevant_gene_pct, max_workers, random_state)
        compute_gene_redundancy(adata, max_workers, random_state)
        compute_gene_complementarity(adata, max_workers, random_state)
        generate_gene_clusters(adata)
    logger.info("Gene clustering done!")


def compute_gene_closeness(adata: ad.AnnData, centers: np.ndarray) -> np.ndarray:
    """
    This function is only used in GeneClust-fast. Firstly, it computes the distance of each gene to its cluster mean,
    min-max normalized in cluster. The gene closeness is then computed as 1 - normalized distances.

    Parameters
    ----------
    adata : AnnData
        The annotated data matrix of shape `n_obs` × `n_vars`.
        Rows correspond to cells and columns to genes.
    centers : ndarray
        Cluster centers of shape `n_var_clusters` × `n_components`.
    Returns
    -------
    all_distances : ndarray
        distances of all genes to their cluster centers.
    """
    all_distances = paired_distances(adata.varm['X_pca'], centers[adata.var['cluster']])
    for gene_cluster in np.unique(adata.var['cluster']):
        gene_cluster_mask = adata.var['cluster'] == gene_cluster
        all_distances[gene_cluster_mask] = 1 - minmax_scale(all_distances[gene_cluster_mask])
    return all_distances


def generate_gene_clusters(adata: ad.AnnData):
    """
    This function is only used in GeneClust-ps. It adopts the algorithm in HDBSCAN to generate gene clusters from an MST,
    and computes the outlier score of each gene based on the GLOSH algorithm.

    Parameters
    ----------
    adata : AnnData
        The annotated data matrix of shape `n_obs` × `n_vars`.
        Rows correspond to cells and columns to genes.
    """
    g1_idx, g2_idx = adata.uns['mst_edges'][:, 0], adata.uns['mst_edges'][:, 1]
    per_gene_relevance = adata.var['relevance'].values
    edge_min_relevance = np.minimum(per_gene_relevance[g1_idx], per_gene_relevance[g2_idx])
    edge_redundancy = adata.varp['redundancy'][g1_idx, g2_idx]
    edge_scales = np.maximum(edge_min_relevance, adata.uns['mst_edges_complm']) / edge_redundancy

    MST = np.hstack((adata.uns['mst_edges'], edge_scales.reshape(-1, 1)))
    MST = MST[np.argsort(MST.T[2]), :]
    single_linkage_tree = label(MST)
    condensed_tree = condense_tree(single_linkage_tree, 3)
    stability_dict = compute_stability(condensed_tree)
    labels, probabilities, stabilities = get_clusters(condensed_tree, stability_dict, "eom", False, False, 0., 0)
    adata.var['outlier_score'], adata.var['cluster'] = outlier_scores(condensed_tree), labels
