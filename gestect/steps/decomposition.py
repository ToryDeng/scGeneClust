# -*- coding: utf-8 -*-
# @Time : 2022/5/10 11:34
# @Author : Tory Deng
# @File : decomposition.py
# @Software: PyCharm
import torch
import anndata as ad
from typing import Literal, Optional
from sklearn.decomposition import TruncatedSVD
from torchnmf.nmf import NMF


def two_way_embedding(adata: ad.AnnData,
                      layer: Optional[str],
                      dr_method: Literal['svd', 'nmf'],
                      n_comps: int = 50) -> None:
    if dr_method == 'svd':
        expr = adata.layers[layer] if layer is not None else adata.X
        model = TruncatedSVD(n_components=n_comps, random_state=2022)
        adata.obsm[dr_method], adata.varm[dr_method] = model.fit_transform(expr), model.components_.T
    elif dr_method == 'nmf':
        expr = torch.from_numpy(adata.layers[layer] if layer is not None else adata.raw.X).cuda()
        model = NMF(expr.shape, rank=n_comps).cuda()
        model.fit(expr)
        adata.obsm[dr_method], adata.varm[dr_method] = model.H.cpu().detach().numpy(), model.W.cpu().detach().numpy()
    else:
        raise NotImplementedError(f"{dr_method} has not been implemented!")
    adata.uns['dr_method'] = dr_method

