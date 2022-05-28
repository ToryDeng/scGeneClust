# -*- coding: utf-8 -*-
# @Time : 2022/5/21 16:12
# @Author : Tory Deng
# @File : utils.py
# @Software: PyCharm
import sys
from typing import Literal

import anndata as ad
import scanpy as sc
from loguru import logger


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


def set_logger(verbose: Literal[0, 1, 2] = 1):
    """
    set the verbosity level.

    :param verbose: integer in [0, 1, 2]. 0: only print warnings and errors; 1: also print info; 2: also print debug message
    :return: None
    """
    assert verbose in (0, 1, 2), ValueError(f"`verbose` can only be 0, 1 or2.")

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
    logger.add(sys.stderr, level=level_dict[verbose], format=formatter)
