import numpy as np
import pandas as pd
import scanpy as sc
import anndata as ad
import steps
from typing import Literal


def select(adata: ad.AnnData,  # row: cells, col: genes
           n_selected_genes: int,
           dr_method: Literal['pca', 'glm-pca', 'umap'],
           n_comps: int,
           similarity: Literal['pearson', 'spearman', 'kendall'],
           n_clusters: int
           ):
    genes_as_samples = adata.to_df().T
    reduced = steps.reduce_dimension(genes_as_samples, dr_method, n_comps)
    dist_mtx = steps.compute_gene_similarity(reduced.T, similarity)
    cluster_labels = steps.clustering_genes(dist_mtx, n_clusters)
    per_cluster_gene_scores = [steps.in_cluster_score(genes_as_samples.loc[cluster_labels == i, :]) for i in range(n_clusters)]
    per_cluster_selected_genes = []
    for i, cluster_gene_scores in enumerate(per_cluster_gene_scores.copy()):
        if cluster_gene_scores.shape[0] >= n_selected_genes // n_clusters:
            top_selected_genes = cluster_gene_scores.iloc[:(n_selected_genes // n_clusters), 0]
            per_cluster_selected_genes.append(top_selected_genes.values)
            per_cluster_gene_scores[i].drop(top_selected_genes.index, inplace=True)
        else:
            per_cluster_selected_genes.append(cluster_gene_scores['index'].values)
            per_cluster_gene_scores[i].drop(cluster_gene_scores['index'].index, inplace=True)
    selected_genes = np.concatenate(per_cluster_selected_genes)
    if selected_genes.shape[0] < n_selected_genes:
        idxs = np.random.choice(len(per_cluster_gene_scores), size=n_selected_genes - selected_genes.shape[0], replace=False)
        remains = [per_cluster_gene_scores[idx].iloc[0, 0] for idx in idxs]
        selected_genes = np.concatenate((selected_genes, remains))
    gene_mask = np.isin(adata.var_names, selected_genes)
    if adata.raw is not  None:
        adata.raw = adata.raw[:, gene_mask].to_adata()
    adata = adata[:, gene_mask]
    return adata


# 1000 cells, 5000 genes after log-norm
test_adata = sc.datasets.blobs(n_variables=5000, n_centers=5, n_observations=1000)
selected = select(test_adata, n_selected_genes=2000, dr_method='pca', n_comps=50, similarity='pearson', n_clusters=200)
print(selected)
