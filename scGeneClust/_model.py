# -*- coding: utf-8 -*-
# @Time : 2022/6/3 16:27
# @Author : Tory Deng
# @File : _model.py
# @Software: PyCharm
from typing import Literal, Optional

import anndata as ad
import numpy as np
from sklearn.cluster import MiniBatchKMeans

import scGeneClust.pp as pp
import scGeneClust.tl as tl
from loguru import logger
from .utils import _check_raw_counts, set_logger, select_from_clusters, prepare_GO


def scGeneClust(
        raw_adata: ad.AnnData,
        mode: Literal['fast', 'hc'] = 'fast',
        n_components: int = 50,
        n_gene_clusters: Optional[int] = None,
        n_cell_clusters: Optional[int] = None,
        verbosity: Literal[0, 1, 2] = 1,
        rlv_threshold: float = 0.01,
        scale: int = 2000,
        random_stat: Optional[int] = None
):
    set_logger(verbosity)
    _check_raw_counts(raw_adata)

    copied = raw_adata.copy()
    # preprocessing
    pp.preprocess(copied, mode)
    pp.reduce_dimension(copied, mode, n_components, random_stat)

    if mode == 'fast':
        km = MiniBatchKMeans(n_clusters=n_gene_clusters, random_state=random_stat)
        copied.var['cluster'] = km.fit_predict(copied.varm['pca'])  # gene clustering
        copied.var['score'] = tl.compute_gene_closeness(copied, km.cluster_centers_)
        tl.handle_single_gene_cluster(copied, mode, random_stat)
        tl.filter_constant_genes(copied)
    else:
        tl.find_high_confidence_cells(copied, n_cell_clusters=n_cell_clusters, random_stat=random_stat)
        tl.filter_low_confidence_cells(copied)
        tl.filter_irrelevant_gene(copied, rlv_threshold, random_stat)
        tl.clustering(copied, scale, random_stat)

    selected_genes = select_from_clusters(copied, mode)

    # TODO: remove this preparation for GO analysis
    # prepare_GO(copied, save='cache/')
    # if 'disease_related' in copied.var:
    #     copied.var.groupby('cluster')['disease_related'].value_counts().to_excel('cache/disease_analysis.xlsx')
    #     logger.info("disease analysis done!")

    # check if all selected features are in var_names
    is_selected = np.isin(raw_adata.var_names, selected_genes)
    if is_selected.sum() != selected_genes.shape[0]:
        raise RuntimeError(
            f"Only found {is_selected.sum()} selected genes in adata.var_names, not {selected_genes.shape[0]}.")
    return selected_genes
