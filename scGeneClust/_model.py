# -*- coding: utf-8 -*-
# @Time : 2022/6/3 16:27
# @Author : Tory Deng
# @File : _model.py
# @Software: PyCharm
import anndata as ad
import numpy as np
import scGeneClust.pp as pp
import scGeneClust.tl as tl
from typing import Literal, Optional
from scGeneClust._utils import _check_raw_counts, _check_params, set_logger, select_from_clusters


def scGeneClust(
        adata: ad.AnnData,
        mode: Literal['one-way', 'two-way'] = 'two-way',
        n_components: int = 50,
        n_gene_clusters: Optional[int] = None,
        n_gene_neighbors: Optional[int] = None,
        n_cell_clusters: Optional[int] = None,
        subset: bool = False,
        verbosity: Literal[0, 1, 2] = 1,
        random_stat: Optional[int] = None,
        **kwargs
):
    set_logger(verbosity)
    _check_raw_counts(adata)
    n_gene_clusters, n_gene_neighbors, n_cell_clusters = \
        _check_params(mode, n_components, n_gene_clusters, n_gene_neighbors, n_cell_clusters, verbosity)

    pp.preprocess(adata, mode)
    pp.reduce_dimension(adata, mode, n_components, random_stat)

    tl.do_clustering(adata, mode, n_gene_clusters, n_gene_neighbors, n_cell_clusters, random_stat, **kwargs)
    filtered_adata = tl.filter(adata, mode)
    if mode == 'two-way':
        tl.score_discriminative_gene(filtered_adata)
        adata.var['score'] = filtered_adata.var['score']
    selected_genes = select_from_clusters(filtered_adata, mode)
    # check if all selected features are in var_names
    is_selected = np.isin(adata.var_names, selected_genes)
    if is_selected.sum() != selected_genes.shape[0]:
        raise RuntimeError(f"Only found {is_selected.sum()} selected genes in adata.var_names, not {selected_genes.shape[0]}.")
    # subset adata or return selected genes
    if subset:
        adata._inplace_subset_var(is_selected)
    else:
        return selected_genes
