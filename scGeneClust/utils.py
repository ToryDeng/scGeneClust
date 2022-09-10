# -*- coding: utf-8 -*-
# @Time : 2022/6/3 16:46
# @Author : Tory Deng
# @File : utils.py
# @Software: PyCharm
import sys
from functools import partial
from itertools import combinations
from multiprocessing import cpu_count
from multiprocessing.pool import Pool
from typing import Literal

import anndata as ad
import numpy as np
import pandas as pd
import scanpy as sc
from loguru import logger
from scipy.spatial.distance import squareform
from sklearn.feature_selection import mutual_info_classif
from sklearn.metrics import pairwise_distances


def _check_raw_counts(adata: ad.AnnData):
    """Check whether the data in `adata.X` is raw counts."""
    if adata.raw is not None:
        raise ValueError("Expect count data in `.raw` as input.")
    else:
        if not np.all(adata.X % 1 == 0):
            raise ValueError("The input data in `.X` may have been normalized.")


def _check_params(raw_adata, version, n_gene_clusters, n_cell_clusters, top_percent_rel, scale, verbosity, random_stat):
    """Check input parameters of the main function."""
    _check_raw_counts(raw_adata)
    if version == 'fast':
        if n_gene_clusters is None:
            raise ValueError("You must specify `n_gene_clusters` in GeneClust-fast.")
        if n_cell_clusters is not None:
            raise UserWarning("The parameter `n_cell_clusters` doesn't need to be specified and is ignored.")
        if top_percent_rel is not None:
            raise UserWarning("The parameter `top_percent_relevance` doesn't need to be specified and is ignored.")
        if scale is not None:
            raise UserWarning("The parameter `scale` doesn't need to be specified and is ignored.")
    elif version == 'ps':
        if n_gene_clusters is not None:
            raise UserWarning("The parameter `n_gene_clusters` doesn't need to be specified and is ignored.")
        if n_cell_clusters is None:
            raise ValueError("You must specify `n_cell_clusters` in GeneClust-ps.")
        if top_percent_rel is None:
            raise ValueError("You must specify `top_percent_relevance` in GeneClust-ps.")
        if scale is None:
            raise ValueError("You must specify `scale` in GeneClust-ps.")
    else:
        raise ValueError("The parameter `version` can only be `fast` or `ps`.")

    if verbosity not in (0, 1, 2):
        raise ValueError("The parameter `verbosity` can only be 0, 1, 2.")
    if random_stat is not None and not isinstance(random_stat, int):
        raise ValueError("The parameter `random_stat` must be None or an integer.")


def _check_all_selected(selected_genes, raw_adata):
    """Check whether all selected genes are in the input `AnnData` object."""
    is_selected = np.isin(raw_adata.var_names, selected_genes)
    if is_selected.sum() != selected_genes.shape[0]:
        msg = f"Only found {is_selected.sum()} selected genes in `adata.var_names`, not {selected_genes.shape[0]}."
        raise RuntimeError(msg)


def load_PBMC3k(min_genes: int = 200, min_cells: int = 3) -> ad.AnnData:
    """
    Load the PBMC3k dataset as an example.

    :param min_genes: Minimum number of genes expressed required for a cell to pass filtering.
    :param min_cells: Minimum number of cells expressed required for a gene to pass filtering.
    :return: The PBMC3k dataset as an `AnnData` object.
    """
    example = sc.datasets.pbmc3k()
    sc.pp.filter_cells(example, min_genes=min_genes)
    sc.pp.filter_genes(example, min_cells=min_cells)
    example.X = example.X.toarray()
    return example


def set_logger(verbosity: Literal[0, 1, 2] = 1):
    """
    Set the verbosity level.

    :param verbosity: 0 (only print warnings and errors), 1 (also print info), 2 (also print debug messages)
    """
    def formatter(record: dict):
        if record['level'].name in ('DEBUG', 'INFO'):
            return "<level>{level: <5}</level> | " \
                   "<level>{message}\n</level>"
        else:
            return "<green>{time:YYYY-MM-DD HH:mm:ss.SSS}</green> | " \
                   "<level>{level: <8}</level> | " \
                   "<cyan>{name}</cyan>:<cyan>{function}</cyan>:<cyan>{line}</cyan> - <level>{message}\n</level>"

    level_dict = {0: 'WARNING', 1: 'INFO', 2: 'DEBUG'}
    logger.remove()
    logger.add(sys.stderr, level=level_dict[verbosity], format=formatter)


def select_from_clusters(
        adata: ad.AnnData,
        version: str,
) -> np.ndarray:
    """
    Select features from gene clusters.

    :param adata: The annotated matrix.
    :param version: The version of GeneClust.
    :return: An ndarray of selected features.
    """
    assert 'cluster' in adata.var, KeyError(f"Column not found in `.var`: 'cluster'")
    assert 'score' in adata.var, KeyError(f"Column not found in `.var`: 'score'")
    df = adata.var.loc[:, ('cluster', 'score')].rename_axis('gene')
    grouped = df.groupby(by='cluster')['score']
    if version == 'fast':
        max_genes = grouped.nlargest(1).reset_index(level=1)['gene'].values
        min_genes = grouped.nsmallest(1).reset_index(level=1)['gene'].values
        selected_features = np.unique(np.concatenate([max_genes, min_genes]))  # the max and min values may be the same
    else:
        selected_features = adata.var_names[adata.var['representative']]
    logger.debug(f"Selected {selected_features.shape[0]} features")
    return selected_features
