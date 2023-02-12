# -*- coding: utf-8 -*-
# @Time : 2022/6/3 16:46
# @Author : Tory Deng
# @File : _utils.py
# @Software: PyCharm
import sys
from typing import Literal, Tuple

import anndata as ad
import cv2
import numpy as np
import scanpy as sc
import squidpy as sq
from loguru import logger


def load_PBMC3k(min_genes: int = 200, min_cells: int = 3) -> ad.AnnData:
    """s
    Load the PBMC3k dataset as an example.

    Parameters
    ----------
    min_genes: int
        Minimum number of genes expressed required for a cell to pass filtering.
    min_cells: int
        Minimum number of cells expressed required for a gene to pass filtering.

    Returns
    -------
    adata : ad.AnnData
        The PBMC3k dataset as an `AnnData` object.
    """
    adata = sc.datasets.pbmc3k()
    sc.pp.filter_cells(adata, min_genes=min_genes)
    sc.pp.filter_genes(adata, min_cells=min_cells)
    adata.X = adata.X.toarray()
    return adata


def load_simulated_data(n_genes: int = 15000, n_celltype: int = 5, n_observations: int = 1000) -> ad.AnnData:
    """
    Gaussian Blobs.

    Parameters
    ----------
    n_genes
        The number of genes.
    n_celltype
        The number of cell types.
    n_observations
        The number of cells.
    Returns
    -------
    adata
        A simulated `AnnData` object.
    """
    adata = sc.datasets.blobs(n_variables=n_genes, n_centers=n_celltype, n_observations=n_observations)
    adata.X[adata.X < 0] = 0
    adata.X[adata.X > 5] = 0
    adata.X = np.round(adata.X, decimals=0)
    adata.obs.rename(columns={'blobs': 'celltype'}, inplace=True)
    sc.pp.filter_cells(adata, min_genes=200)
    sc.pp.filter_genes(adata, min_cells=3)
    return adata


def load_mouse_brain(min_genes: int = 200, min_spots: int = 3) -> Tuple[ad.AnnData, np.ndarray]:
    adata = sq.datasets.visium('V1_Adult_Mouse_Brain', include_hires_tiff=True)
    adata.var_names_make_unique()
    img = cv2.imread(adata.uns['spatial']['V1_Adult_Mouse_Brain']['metadata']['source_image_path'])
    sc.pp.filter_cells(adata, min_genes=min_genes)
    sc.pp.filter_genes(adata, min_cells=min_spots)
    adata.X = adata.X.toarray()
    return adata, img


def set_logger(verbosity: Literal[0, 1, 2] = 1):
    """
    Set the verbosity level.

    Parameters
    ----------
    verbosity
        0 (only print warnings and errors), 1 (also print info), 2 (also print debug messages)
    """
    def formatter(record: dict):
        if record['level'].name in ('DEBUG', 'INFO'):
            return "<green>{time:YYYY-MM-DD HH:mm:ss.SSS}</green> | " \
                   "<level>{level: <5}</level> | " \
                   "<level>{message}\n</level>"
        else:
            return "<green>{time:YYYY-MM-DD HH:mm:ss.SSS}</green> | " \
                   "<level>{level: <8}</level> | " \
                   "<cyan>{name}</cyan>:<cyan>{function}</cyan>:<cyan>{line}</cyan> - <level>{message}\n</level>"

    level_dict = {0: 'WARNING', 1: 'INFO', 2: 'DEBUG'}
    logger.remove()
    logger.add(sys.stdout, colorize=True, level=level_dict[verbosity], format=formatter)
