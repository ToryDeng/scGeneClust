# -*- coding: utf-8 -*-
# @Time : 2022/5/10 11:52
# @Author : Tory Deng
# @File : functions.py
# @Software: PyCharm
import anndata as ad
from typing import Literal, Optional
from queue import Queue
from threading import Thread
import gestect.steps as steps
from GeneClust.utils import subset_adata
from GeneClust.functions import iteratively_select


def two_way_clustering(adata: ad.AnnData,
                       n_cell_clusters: int,
                       gene_clustering_method: str,
                       n_gene_clusters: int,
                       confidence: Literal['single_proba', 'mnn_proba', 'mnn_graph']):
    q = Queue()
    t1 = Thread(target=steps.find_confident_cells, args=(adata.obsm[adata.uns['dr_method']], n_cell_clusters, q, confidence))
    t2 = Thread(target=steps.cluster_genes, args=(adata.varm[adata.uns['dr_method']], gene_clustering_method, n_gene_clusters, q))

    t1.start()
    t2.start()
    t1.join()
    t2.join()

    cell_result = q.get()
    adata.obs['cluster'], adata.obs['confident'] = cell_result['cluster'], cell_result['confident']
    print(f"removing {adata.n_obs - adata.obs['confident'].sum()} unconfident cells...")
    adata = adata[adata.obs['confident'], :]

    adata.var['cluster'] = q.get()
    top_large = adata.var['cluster'].value_counts(normalize=True).head(5)
    print("Percentage of the top 5 largest gene clusters:", top_large, sep='\n')

    return adata


def gestect(adata:ad.AnnData,
            n_selected_genes: int,
            use_rep: Optional[str] = None,
            dr_method: Literal['svd', 'nmf'] = 'svd',
            n_components: int = 50,
            confidence: Literal['single_proba', 'mnn_proba', 'mnn_graph'] = 'single_proba',
            n_cell_clusters: int = 2,
            gene_clustering: Literal['ms', 'agg', 'gmm'] = 'gmm',
            n_gene_clusters: int = 300,
            gene_score: Literal['f_stat'] = 'f_stat',
            top_n_genes: int = 1,
            gene_cluster_score: Literal['top_mean'] = 'top_mean',
            return_genes: bool = False
            ):
    steps.two_way_embedding(adata, use_rep, dr_method, n_components)
    adata = two_way_clustering(adata, n_cell_clusters, gene_clustering, n_gene_clusters, confidence)
    steps.compute_gene_score(adata, gene_score)
    steps.compute_gene_cluster_score(adata, top_n_genes, gene_cluster_score)
    adata.var.to_csv(f"cache/{adata.uns['data_name']}_var.csv")
    selected_genes = iteratively_select(adata, n_selected_genes, 'cluster', 'significant', 'gene_score')
    # only preserve selected genes in adata
    filtered_adata = subset_adata(adata, selected_genes, inplace=False)
    return filtered_adata.var_names.to_frame(index=False, name='Gene') if return_genes else filtered_adata

