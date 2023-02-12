# -*- coding: utf-8 -*-
# @Time : 2023/1/5 1:54
# @Author : Tory Deng
# @File : _validation.py
# @Software: PyCharm
import os
from typing import Literal, Optional

import anndata as ad
import numpy as np
from loguru import logger
from scipy.sparse import issparse


def check_args(
        adata: ad.AnnData,
        img: np.ndarray,
        version: Literal['fast', 'ps'],

        n_gene_clusters: int,
        n_cell_clusters: int,
        n_components: int,
        relevant_gene_pct: int,
        post_hoc_filtering: bool,
        modality: Literal['sc', 'st'],
        shape: Literal['hexagon', 'square'],
        return_info: bool,
        subset: bool,
        max_workers: Optional[int],
        verbosity: Literal[0, 1, 2],
        random_state: int,
):
    """Check all parameters"""
    # check base args
    if verbosity not in (0, 1, 2):
        raise ValueError(f"Expected `verbosity` in {{0, 1, 2}}, got {verbosity}.")
    if version not in ('fast', 'ps'):
        raise ValueError(f"Expected `version` in {{'fast', 'ps'}}, got {version}.")
    if modality not in ('sc', 'st'):
        raise ValueError(f"Expected `modality` in {{'sc', 'st'}}, got {modality}.")

    if not isinstance(random_state, int):
        raise TypeError(f"Expected `random_stat` to be an integer, got {type(random_state)}.")

    if not isinstance(max_workers, int):
        raise TypeError(f"Expected `max_workers` to be a positive integer, got {type(max_workers)}.")
    if max_workers < 1:
        raise ValueError(f"Expected `max_workers` to be a positive integer, got {max_workers}.")
    if max_workers > os.cpu_count():
        raise ValueError(f"Worker limit exceeded. Maximum {os.cpu_count()}, got {max_workers}.")

    if not isinstance(return_info, bool):
        raise TypeError(f"Expected `return_info` to be bool, got {type(return_info)}.")
    if not isinstance(subset, bool):
        raise TypeError(f"Expected `subset` to be bool, got {type(subset)}.")

    check_raw_counts(adata)
    check_gene_clustering(img, version, n_gene_clusters, n_cell_clusters, n_components, relevant_gene_pct,
                          post_hoc_filtering, modality, shape)


def check_raw_counts(adata: ad.AnnData):
    """Check whether `adata.X` contains raw counts"""
    if isinstance(adata, ad.AnnData):
        if adata.raw is not None:
            raise ValueError("Expected raw counts in `adata.raw` as input.")
        else:
            X = adata.X.toarray() if issparse(adata.X) else adata.X
            if not np.all(X % 1 == 0):
                raise ValueError("Expected raw counts, got normalized counts possibly.")
    else:
        raise TypeError(f"Expected an `AnnData` object, got {type(adata)}.")


def check_gene_clustering(
        img: np.ndarray,
        version: Literal['fast', 'ps'],
        n_gene_clusters: int,
        n_obs_clusters: int,
        n_components: int,
        relevant_gene_pct: int,
        post_hoc_filtering: bool,
        modality: Literal['sc', 'st'],
        shape: Literal['hexagon', 'square'],
):
    """Check params used in gene clustering"""
    # check post_hoc_filtering
    if not isinstance(post_hoc_filtering, bool):
        raise TypeError(f"Expected `post_hoc_filtering` to be bool, got {type(post_hoc_filtering)}.")

    if version == 'fast':  # GeneClust-fast requires `n_clusters` to be given
        if not isinstance(n_gene_clusters, int):
            raise TypeError(f"Expected `n_gene_clusters` to be an integer, got {type(n_gene_clusters)}.")
        else:
            if n_gene_clusters <= 1:
                raise ValueError(f"Expected `n_gene_clusters` > 1, got {n_gene_clusters}.")

        if modality == 'st':
            raise ValueError(f"GeneClust-fast does not support spatial transcriptomics. Please set `version == 'ps'`.")
    else:
        if n_gene_clusters is not None:  # GeneClust-ps does not require `n_var_clusters` to be given
            raise TypeError(f"Expected `n_gene_clusters` to be None, got {type(n_gene_clusters)}.")
        if modality == 'st':
            if img is None:
                raise RuntimeWarning(f"Image is not available. This could lower the accuracy of spaGCN.")
            elif not isinstance(img, np.ndarray):
                raise TypeError(f"Expected `image` to be an ndarray, got {type(img)}.")
            else:
                pass
        if shape not in ('hexagon', 'square'):
            raise ValueError(f"Expected `shape` in {{'hexagon', 'square'}}, got {shape}.")
        # check n_obs_clusters
        if not isinstance(n_obs_clusters, int):
            raise TypeError(f"Expected `n_obs_clusters` to be an integer, got {type(n_obs_clusters)}.")
        else:
            if n_obs_clusters <= 1:
                raise ValueError(f"Expected `n_obs_clusters` > 1, got {n_obs_clusters}.")
        # check n_components
        if not isinstance(n_components, int):
            raise TypeError(f"Expected `n_components` to be an integer, got {type(n_components)}.")
        else:
            if n_components <= 1:
                raise ValueError(f"Expected `n_components` > 1, got {n_components}.")
        # check relevant_gene_pct
        if not isinstance(relevant_gene_pct, int):
            raise TypeError(f"Expected `relevant_gene_pct` to be an integer, got {type(relevant_gene_pct)}.")
        else:
            if relevant_gene_pct <= 0 or relevant_gene_pct >= 100:
                raise ValueError(f"Expected `relevant_gene_pct` between (0, 100), got {relevant_gene_pct}.")


def check_all_genes_selected(adata: ad.AnnData, selected_genes: np.ndarray):
    """Check whether all selected genes are in the input `AnnData` object."""
    is_selected = np.isin(adata.var_names, selected_genes)
    if (n_selected_genes := is_selected.sum()) != selected_genes.shape[0]:
        msg = f"Found only {n_selected_genes} selected genes in `adata.var_names`, not {selected_genes.shape[0]}."
        raise RuntimeError(msg)
    logger.opt(colors=True).info(f"Selected <yellow>{n_selected_genes}</yellow> genes.")


