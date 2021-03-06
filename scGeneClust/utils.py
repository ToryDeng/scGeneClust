# -*- coding: utf-8 -*-
# @Time : 2022/6/3 16:46
# @Author : Tory Deng
# @File : utils.py
# @Software: PyCharm
import os
import sys
from typing import Literal, Optional

import anndata as ad
import numpy as np
import scanpy as sc
from rpy2.robjects import globalenv, r
import anndata2ri
from loguru import logger


def _check_raw_counts(adata: ad.AnnData):
    if adata.raw is not None:
        raise ValueError("Expect count data in `.raw` as input.")
    else:
        if not np.all(adata.X % 1 == 0):
            raise ValueError("The input data in `.X` may have been normalized.")


def load_example_adata(min_genes: int = 200, min_cells: int = 3) -> ad.AnnData:
    example = sc.datasets.pbmc3k()
    sc.pp.filter_cells(example, min_genes=min_genes)
    sc.pp.filter_genes(example, min_cells=min_cells)
    example.X = example.X.toarray()
    example.var['original_gene'] = example.var_names
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
    Set the verbosity level.

    Parameters
    ----------
    verbosity
      integer in [0, 1, 2]. 0: only print warnings and errors; 1: also print info; 2: also print debug message
    Returns
    -------

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
    assert 'cluster' in adata.var, KeyError(f"Column not found in `.var`: 'cluster'")
    assert 'score' in adata.var, KeyError(f"Column not found in `.var`: 'score'")
    df = adata.var.loc[:, ('cluster', 'score')].rename_axis('gene')
    grouped = df.groupby(by='cluster')['score']
    if mode == 'fast':
        max_genes = grouped.nlargest(1).reset_index(level=1)['gene'].values
        min_genes = grouped.nsmallest(1).reset_index(level=1)['gene'].values
        selected_features = np.unique(np.concatenate([max_genes, min_genes]))  # the max and min values may be the same
    else:
        selected_features = adata.var.index
    logger.debug(f"Selected {selected_features.shape[0]} features")
    return selected_features


def prepare_GO(raw_adata: ad.AnnData, save: Optional[str] = None, name=None, save_type='Rdata'):
    assert 'original_gene' in raw_adata.var, ValueError("Column 'original_gene' not in adata.var")
    assert 'cluster' in raw_adata.var and 'score' in raw_adata.var, ValueError("Must run `scGeneClust` first!")

    sc.pp.highly_variable_genes(raw_adata, n_top_genes=500, flavor='seurat_v3')
    if save is not None:
        if not os.path.exists(save):
            os.makedirs(save)

        if name is None:
            if 'data_name' not in raw_adata.uns:
                name = 'data'
            else:
                name = raw_adata.uns['data_name']
        if save_type == 'csv':
            path_to_save = os.path.join(save, f"{name}.csv")
            raw_adata.var.loc[:, ['original_gene', 'cluster', 'score', 'variances_norm']].to_csv(path_to_save,
                                                                                                 index=False)
        else:
            path_to_save = os.path.join(save, f"{name}.Rdata")
            anndata2ri.activate()
            globalenv['sce'], globalenv['path'] = anndata2ri.py2rpy(raw_adata), path_to_save
            r("""
            save(sce, file = path)
            """)
            anndata2ri.deactivate()

        logger.info(f"data has been saved at {path_to_save}")
    logger.info("Preparation done!")
