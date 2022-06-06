# -*- coding: utf-8 -*-
# @Time : 2022/6/3 16:46
# @Author : Tory Deng
# @File : _utils.py
# @Software: PyCharm
import anndata as ad
import numpy as np
from loguru import logger
import scanpy as sc
import sys
from numpy.random import default_rng
from typing import Literal, Optional


def _check_raw_counts(adata: ad.AnnData):
    if adata.raw is not None:
        raise ValueError("Expect count data in `.raw` as input.")
    else:
        if not np.all(adata.X % 1 == 0):
            raise ValueError("The input data in `.X` seems to have been normalized.")


def _check_params(
        mode: Literal['one-way', 'two-way'],
        n_components: int,
        n_gene_clusters: Optional[int],
        n_gene_neighbors: Optional[int],
        n_cell_clusters: Optional[int],
        verbosity: Literal[0, 1, 2]
):
    assert isinstance(n_components, int) and n_components > 0, TypeError("`n_components` must be a nonnegative integer!")
    assert verbosity in (0, 1, 2), ValueError(f"`verbose` can only be 0, 1 or 2.")

    if mode == 'one-way':
        if n_gene_clusters is None:
            logger.info("`n_gene_clusters` is None. Using default value (200)")
            n_gene_clusters = 200
        logger.info("Will ignore `n_gene_neighbors`")
        logger.info("Will ignore `n_cell_clusters`")
        n_gene_neighbors, n_cell_clusters = None, None
    elif mode == 'two-way':
        logger.info("Will ignore `n_gene_clusters`")
        n_gene_clusters = None
        if n_gene_neighbors is None:
            logger.info("`n_gene_neighbors` is None. Using default value (30)")
            n_gene_neighbors = 30
        assert n_cell_clusters is not None, \
            ValueError("Argument `n_cell_clusters` must be specified in 'two-way' mode!")
    else:
        raise ValueError(f"Argument `mode` can only be 'one-way' or 'two-way', not '{mode}'.")
    return n_gene_clusters, n_gene_neighbors, n_cell_clusters


def load_example_adata(min_genes: int = 200, min_cells: int = 3) -> ad.AnnData:
    example = sc.datasets.pbmc3k()
    sc.pp.filter_cells(example, min_genes=min_genes)
    sc.pp.filter_genes(example, min_cells=min_cells)
    example.X = example.X.toarray()
    example.raw = example
    sc.pp.normalize_total(example)
    example.layers['normalized'] = example.X.copy()
    sc.pp.log1p(example)
    example.layers['log-normalized'] = example.X.copy()
    sc.pp.scale(example)
    example.uns['data_name'] = 'pbmc3k'
    return example


def set_logger(verbosity: Literal[0, 1, 2] = 1):
    """
    set the verbosity level.

    :param verbosity: integer in [0, 1, 2]. 0: only print warnings and errors; 1: also print info; 2: also print debug message
    :return: None
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
        mode: str,
) -> np.ndarray:
    """
    Select 2 genes from each gene clusters.

    :param adata: The AnnData object. adata.var must contain 'cluster' and 'score' columns and is indexed by genes
    :param mode:
    :return: ndarray, selected features
    """
    assert 'cluster' in adata.var, KeyError(f"Column not found in `.var`: 'cluster'")
    assert 'score' in adata.var, KeyError(f"Column not found in `.var`: 'score'")
    df = adata.var.loc[:, ('cluster', 'score')].rename_axis('gene')
    grouped = df.groupby(by='cluster')['score']
    if mode == 'one-way':
        max_genes = grouped.nlargest(1).reset_index(level=1)['gene'].values
        min_genes = grouped.nsmallest(1).reset_index(level=1)['gene'].values
        selected_features = np.unique(np.concatenate([min_genes, max_genes]))  # the max and min values may be the same
    else:
        selected_features = grouped.nlargest(2).reset_index(level=1)['gene'].values
    logger.debug(f"Selected {selected_features.shape[0]} features")
    return selected_features


