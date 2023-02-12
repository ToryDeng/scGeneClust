# -*- coding: utf-8 -*-
# @Time : 2022/6/3 16:27
# @Author : Tory Deng
# @File : _model.py
# @Software: PyCharm
import os
from typing import Literal, Optional, Union, Tuple

import anndata as ad
import numpy as np
from loguru import logger
from scipy.sparse import issparse

import scGeneClust.pp as pp
import scGeneClust.tl as tl
from ._utils import set_logger
from ._validation import check_args, check_all_genes_selected


def scGeneClust(
        raw_adata: ad.AnnData,
        image: np.ndarray = None,
        n_var_clusters: int = None,
        n_obs_clusters: int = None,
        n_components: int = 10,
        relevant_gene_pct: int = 20,
        post_hoc_filtering: bool = True,
        version: Literal['fast', 'ps'] = 'fast',
        modality: Literal['sc', 'st'] = 'sc',
        shape: Literal['hexagon', 'square'] = 'hexagon',
        return_info: bool = False,
        subset: bool = False,
        max_workers: int = os.cpu_count() - 1,
        verbosity: Literal[0, 1, 2] = 1,
        random_state: int = 0
) -> Optional[Union[Tuple[ad.AnnData, np.ndarray], np.ndarray]]:
    """
    This function is the common interface for *GeneClust-fast* and *GeneClust-ps*.

    Parameters
    ----------
    raw_adata : AnnData
        The annotated data matrix of shape `n_obs` × `n_vars`.
        Rows correspond to cells and columns to genes.
        The raw counts should be in `raw_adata.X`.
    image : ndarray
        The image of tissue section.
    n_var_clusters : int
        The number of clusters in gene clustering. Only valid in GeneClust-fast.
    n_obs_clusters : int
        The number of clusters in cell clustering used to find high-confidence cells. Only valid in GeneClust-ps.
    n_components : int, default=10
        The number of principal components used along with the first component. Only valid in GeneClust-ps.
    relevant_gene_pct: int, default=20
        The percentage of relevant genes. This parameter should be between 0 and 100. Only valid in GeneClust-ps.
    post_hoc_filtering : bool, default=True
        Whether to find outliers in singleton gene clusters (in GeneClust-fast) or low-density genes (in GeneClust-ps)
        after gene clustering.
    version : Literal['fast', 'ps'], default='fast'
        Choose the version of GeneClust.
    modality : Literal['sc', 'st'], default='sc'
        Type of the dataset. 'sc' for scRNA-seq data, 'st' for spatially resolved transcriptomics (SRT) data.
    shape : Literal['hexagon', 'square'], default='hexagon'
        The shape of spot neighbors. 'hexagon' for Visium data, 'square' for ST data.
    return_info: bool, default=False
        If `False`, only return names of selected genes.
        Otherwise, return an `AnnData` object which contains intermediate results generated during feature selection.
    subset: bool, default=False
        If `True`, inplace subset to selected genes otherwise merely return the names of selected genes
        (and intermediate results recorded in an `AnnData` object, depending on the value of `return_info`).
    max_workers : int, default=os.cpu_count() - 1
        The maximum value of workers which can be used during feature selection. Default is the number of CPUs - 1.
    verbosity : Literal[0, 1, 2], default=1
        The verbosity level.
        If 0, only prints warnings and errors.
        If 1, prints info-level messages, warnings and errors.
        If 2, prints debug-level messages, info-level messages, warnings and errors.
    random_state : int, default=0
        Change to use different initial states for the optimization.

    Returns
    -------
    Depending on `subset` and `return_info`, returns names of selected genes (and intermediate results),
    or inplace subsets to selected genes and returns `None`.

    copied_adata : AnnData
        Stores intermediate results generated during feature selection.
        The normalized counts are stored in `copied_adata.layers['pearson_norm']`.
        The cell-level principal components are stored in `copied_adata.varm['X_pca']`.
        The gene cluster labels are in `copied_adata.var['cluster']`.
        For GeneClust-fast, the closeness of genes to their cluster centers are in `copied_adata.var['closeness']`.
        For GeneClust-ps, the gene-level principal components are in `copied_adata.obsm['X_pca']`.
        The high-confidence cell cluster labels are in `copied_adata.obs['cluster']`.
        Low-confidence cell clusters are filtered.
        Genes relevance values are in `copied_adata.var['relevance']`. Irrelevant genes are filtered.
        Gene redundancy values are in `copied_adata.varp['redundancy']`.
        MST edges are in `copied_adata.uns['mst_edges']` as an ndarray of shape (n_edges, 2).
        Gene complementarity values are in `copied_adata.uns['mst_edges_complm']` as an ndarray of shape (n_edges, ).
        Representative genes are indicated by `copied_adata.var['representative']`.
    selected_genes : ndarray
        Names of selected genes.s

    Examples
    -------
    >>> from scGeneClust import scGeneClust, load_PBMC3k
    >>>
    >>>
    >>> adata = load_PBMC3k()
    >>> selected_genes_fast = scGeneClust(adata, version='fast', n_var_clusters=200)
    >>> selected_genes_ps = scGeneClust(adata, version='ps', n_obs_clusters=7)
    """
    # check arguments
    check_args(raw_adata, image, version, n_var_clusters, n_obs_clusters, n_components, relevant_gene_pct,
               post_hoc_filtering, modality, shape, return_info, subset, max_workers, verbosity, random_state)
    # set log level
    set_logger(verbosity)
    # feature selection starts
    logger.opt(colors=True).info(
        f"Performing <magenta>GeneClust-{version}</magenta> "
        f"on <magenta>{'scRNA-seq' if modality == 'sc' else 'SRT'}</magenta> data, "
        f"with <yellow>{max_workers}</yellow> workers."
    )
    copied_adata = raw_adata.copy()
    copied_adata.X = raw_adata.X.toarray() if issparse(raw_adata.X) else raw_adata.X

    # preprocessing
    pp.normalize(copied_adata, modality)
    pp.reduce_dim(copied_adata, version, random_state)
    # gene clustering
    tl.cluster_genes(copied_adata, image, version, modality, shape, n_var_clusters, n_obs_clusters, n_components, relevant_gene_pct, max_workers, random_state)
    # select features from gene clusters
    selected_genes = tl.select_from_clusters(copied_adata, version, post_hoc_filtering, random_state)
    check_all_genes_selected(raw_adata, selected_genes)

    if subset:
        raw_adata._inplace_subset_var(selected_genes)
        logger.opt(colors=True).info(f"<magenta>GeneClust-{version}</magenta> finished.")
        return None

    logger.opt(colors=True).info(f"GeneClust-{version} finished.")
    if return_info:
        return copied_adata, selected_genes
    else:
        return selected_genes
