# -*- coding: utf-8 -*-
# @Time : 2022/4/8 21:03
# @Author : Tory Deng
# @File : score.py
# @Software: PyCharm
"""
1. in-cluster scoring (how to rank genes in each cluster)
2. (optional) inter-cluster scoring (how to rank clusters)
"""
import pandas as pd
import anndata2ri
import anndata as ad
import traceback
from typing import Optional, Literal
from rpy2.robjects import r, globalenv
from rpy2.robjects.packages import importr
from GeneClust.utils import HiddenPrints


def in_cluster_score(adata: ad.AnnData, score: Literal['var', 'm3drop']):
    """
    Compute the scores for each gene in a cluster.

    :param adata: The anndata object.
    :param score: Which kind of score to use.
    """
    if score == 'var':
        adata.var[score] = adata.X.var(axis=0)
        # return cluster_expr.var(axis=1).sort_values(ascending=False).reset_index()
    elif score == 'm3drop':
        adata.var[score] = M3Drop_compute_importance(adata)['importance']
    else:
        raise NotImplementedError(f"{score} has not been implemented!")


def M3Drop_compute_importance(adata: ad.AnnData) -> Optional[pd.DataFrame]:
    """
    Select features by M3Drop. The input anndata object contains norm and raw data. The raw data is normalized in R.

    :param adata: anndata object.
    :return: A dataframe. The first column contains gene names, and the second column contains 1 - p.value.
    """
    try:
        with HiddenPrints():
            anndata2ri.activate()
            importr('M3Drop')
            globalenv['sce'] = anndata2ri.py2rpy(adata.raw.to_adata())
            r("""
            norm <- M3DropConvertData(assay(sce, 'X'), is.counts=TRUE)
            DE_genes <- M3DropFeatureSelection(norm, mt_method="fdr", mt_threshold=1, suppress.plot = TRUE)
            """)
            result = r("DE_genes").drop(columns=['effect.size', 'q.value', 'Gene'])
            result['importance'] = 1 - result['p.value']
            anndata2ri.deactivate()
            return result
    except:
        traceback.print_exc()
