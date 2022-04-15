# -*- coding: utf-8 -*-
# @Time : 2022/4/11 17:33
# @Author : Tory Deng
# @File : utils.py
# @Software: PyCharm
import sys

import anndata as ad
from typing import Union
import numpy as np
import os
import pandas as pd


def subset_adata(adata: ad.AnnData, selected_genes: Union[np.ndarray, pd.Index], inplace=False):
    if isinstance(selected_genes, pd.Index):
        selected_genes = selected_genes.to_numpy()
    gene_mask = adata.var_names.isin(selected_genes)
    if inplace:
        if adata.raw is not None:
            adata.raw = adata.raw[:, gene_mask].to_adata()
            if adata.raw.shape[1] != selected_genes.shape[0]:
                raise RuntimeError(f"{adata.raw.shape[1]} genes in raw data were selected, "
                                   f"not {selected_genes.shape[0]} genes. Please check the gene names.")
        adata._inplace_subset_var(selected_genes)
        if adata.shape[1] != selected_genes.shape[0]:
            raise RuntimeError(
                f"{adata.shape[1]} in norm data were selected, "
                f"not {selected_genes.shape[0]} genes. Please check the gene names."
            )
    else:
        copied_adata = adata.copy()
        subset_adata(copied_adata, selected_genes, inplace=True)
        return copied_adata


class HiddenPrints:
    """
    Hide prints from terminal
    """
    def __enter__(self):
        self._original_stdout = sys.stdout
        self._original_stderr = sys.stderr
        sys.stdout = open(os.devnull, 'w')
        sys.stderr = open(os.devnull, 'w')

    def __exit__(self, exc_type, exc_val, exc_tb):
        sys.stdout.close()
        sys.stderr.close()
        sys.stdout = self._original_stdout
        sys.stderr = self._original_stderr