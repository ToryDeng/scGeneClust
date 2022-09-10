# -*- coding: utf-8 -*-
# @Time : 2022/6/3 16:27
# @Author : Tory Deng
# @File : _model.py
# @Software: PyCharm
from typing import Literal, Optional

import anndata as ad

import scGeneClust.pp as pp
import scGeneClust.tl as tl
from .utils import set_logger, select_from_clusters, _check_params, _check_all_selected


def scGeneClust(
        raw_adata: ad.AnnData,
        version: Literal['fast', 'ps'] = 'fast',
        n_gene_clusters: Optional[int] = None,
        n_cell_clusters: Optional[int] = None,
        top_percent_relevance: Optional[int] = None,
        scale: Optional[int] = None,
        verbosity: Literal[0, 1, 2] = 1,
        random_stat: Optional[int] = None
):
    """
    The main function of GeneClust.

    :param raw_adata: The annotated matrix. GeneClust expects raw counts.
    :param version: The version of GeneClust.
    :param n_gene_clusters: The number of gene clusters. Only used in GeneClust-fast.
    :param n_cell_clusters: The number of cell clusters. Only used in GeneClust-ps.
    :param top_percent_relevance: What percentage of genes with top relevance should be preserved.
    :param scale: The scale factor used in the partition of MST.
    :param verbosity: The verbosity level.
    :param random_stat: Change to use different initial states for the optimization.
    :return: An ndarray of selected features.
    """
    set_logger(verbosity)
    _check_params(raw_adata, version, n_gene_clusters, n_cell_clusters, top_percent_relevance, scale,
                  verbosity, random_stat)

    copied = raw_adata.copy()
    copied.raw = raw_adata

    # preprocessing
    pp.preprocess(copied, version)
    pp.reduce_dimension(copied, version, random_stat)

    # gene clustering
    if version == 'fast':
        tl.gene_clustering_mbkmeans(copied, n_gene_clusters, random_stat)
        tl.handle_single_gene_cluster(copied, version, random_stat)
        tl.filter_constant_genes(copied)
    else:
        tl.find_high_confidence_cells(copied, n_cell_clusters, random_stat)
        tl.filter_low_confidence_cells(copied)
        tl.filter_irrelevant_gene(copied, top_percent_relevance, random_stat)
        tl.gene_clustering_graph(copied, scale, random_stat)
        tl.handle_single_gene_cluster(copied, version, random_stat)

    # select features from gene clusters
    selected_genes = select_from_clusters(copied, version)
    _check_all_selected(selected_genes, raw_adata)
    return selected_genes
