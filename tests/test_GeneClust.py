# -*- coding: utf-8 -*-
# @Time : 2022/5/25 12:51
# @Author : Tory Deng
# @File : test_GeneClust.py
# @Software: PyCharm
import anndata as ad
import numpy as np

import scGeneClust as gc


def test_load_example_adata():
    adata = gc.load_PBMC3k()
    assert adata is not None
    assert isinstance(adata, ad.AnnData)


def test_singl_cell_GeneClust():
    raw_adata = gc.load_PBMC3k()
    info, genes_fast = gc.scGeneClust(raw_adata, n_var_clusters=200, version='fast', verbosity=2, return_info=True)
    info, genes_ps = gc.scGeneClust(raw_adata, n_obs_clusters=7, relevant_gene_pct=5, version='ps', verbosity=2, return_info=True)
    assert genes_fast.shape[0] > 0 and genes_ps.shape[0] > 0
    assert np.all(info.varp['redundancy'] == info.varp['redundancy'].T)


def test_spatial_GeneClust():
    adata, img = gc.load_mouse_brain()
    info, genes_ps = gc.scGeneClust(adata, img, n_obs_clusters=5, version='ps', modality='st', relevant_gene_pct=5, verbosity=2, return_info=True)
    assert genes_ps.shape[0] > 0


