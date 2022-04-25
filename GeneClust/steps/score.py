# -*- coding: utf-8 -*-
# @Time : 2022/4/8 21:03
# @Author : Tory Deng
# @File : score.py
# @Software: PyCharm
"""
1. intra-cluster scoring (how to rank genes in each cluster)
2. inter-cluster scoring (how to rank clusters)
"""
import traceback
from typing import Optional, Literal

import anndata as ad
import scanpy as sc
import anndata2ri
import pandas as pd
from rpy2.robjects import r, globalenv
from rpy2.robjects.packages import importr
from sklearn.metrics import silhouette_samples
from GeneClust.utils import HiddenPrints
import numpy as np
from scipy.sparse import *
from scipy.spatial.distance import euclidean
from skfeature.utility.construct_W import construct_W


def in_cluster_score(adata: ad.AnnData, score: Literal['m3drop', 'seurat', 'center']):
    """
    Compute the scores for each gene in a cluster.

    :param adata: The anndata object.
    :param score: Which kind of score to use.
    """
    if score == 'm3drop':
        adata.var[score] = M3Drop_compute_importance(adata)['importance']
    elif score == 'seurat':
        adata.var[score] = sc.pp.highly_variable_genes(
            adata.raw.to_adata(), n_top_genes=adata.n_vars // 2, flavor='seurat_v3', inplace=False
        )['variances_norm']
    elif score == 'center':
        compute_distance = np.vectorize(euclidean, signature='(n),(n)->()')
        adata.var[score] = compute_distance(
            adata.varm[adata.uns['dr_method']], adata.uns['cluster_centroid'][adata.var['cluster_label']]
        )
    else:
        raise NotImplementedError(f"{score} has not been implemented!")
    adata.uns['intra_cluster_score'] = score


def inter_cluster_score(adata: ad.AnnData, inter_score: Literal['top3', 'silhouette']):
    if inter_score == 'silhouette':
        sample_silhouette_values = silhouette_samples(
            X=adata.varp[adata.uns['distance']], labels=adata.var['cluster_label'], metric='precomputed'
        )
        for cluster in adata.var['cluster_label'].unique():
            cluster_mask = adata.var['cluster_label'] == cluster
            adata.var.loc[cluster_mask, 'cluster_score'] = sample_silhouette_values[cluster_mask].mean()
    elif inter_score == 'top3':
        for cluster in adata.var['cluster_label'].unique():
            cluster_mask = adata.var['cluster_label'] == cluster
            sorted_scores = adata.var.loc[cluster_mask, adata.uns['intra_cluster_score']].sort_values(ascending=False)
            if len(sorted_scores) <= 3:
                adata.var.loc[cluster_mask, 'cluster_score'] = sorted_scores.mean()
            else:
                adata.var.loc[cluster_mask, 'cluster_score'] = sorted_scores.iloc[:3].mean()

    adata.var['use_cluster'] = adata.var['cluster_score'] >= adata.var['cluster_score'].quantile(0.1)


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


def fisher_score(X, y):
    """
    This function implements the fisher score feature selection, steps are as follows:
    1. Construct the affinity matrix W in fisher score way
    2. For the r-th feature, we define fr = X(:,r), D = diag(W*ones), ones = [1,...,1]', L = D - W
    3. Let fr_hat = fr - (fr'*D*ones)*ones/(ones'*D*ones)
    4. Fisher score for the r-th feature is score = (fr_hat'*D*fr_hat)/(fr_hat'*L*fr_hat)-1

    Input
    -----
    X: {numpy array}, shape (n_samples, n_features)
        input data
    y: {numpy array}, shape (n_samples,)
        input class labels

    Output
    ------
    score: {numpy array}, shape (n_features,)
        fisher score for each feature

    Reference
    ---------
    He, Xiaofei et al. "Laplacian Score for Feature Selection." NIPS 2005.
    Duda, Richard et al. "Pattern classification." John Wiley & Sons, 2012.
    """

    # Construct weight matrix W in a fisherScore way
    kwargs = {"neighbor_mode": "supervised", "fisher_score": True, "y": y}
    W = construct_W(X, **kwargs)

    # build the diagonal D matrix from affinity matrix W
    D = np.array(W.sum(axis=1))
    L = W
    tmp = np.dot(np.transpose(D), X)
    D = diags(np.transpose(D), [0])
    Xt = np.transpose(X)
    t1 = np.transpose(np.dot(Xt, D.todense()))
    t2 = np.transpose(np.dot(Xt, L.todense()))
    # compute the numerator of Lr
    D_prime = np.sum(np.multiply(t1, X), 0) - np.multiply(tmp, tmp) / D.sum()
    # compute the denominator of Lr
    L_prime = np.sum(np.multiply(t2, X), 0) - np.multiply(tmp, tmp) / D.sum()
    # avoid the denominator of Lr to be 0
    D_prime[D_prime < 1e-12] = 10000
    lap_score = 1 - np.array(np.multiply(L_prime, 1 / D_prime))[0, :]

    # compute fisher score from laplacian score, where fisher_score = 1/lap_score - 1
    score = 1.0 / lap_score - 1

    return np.transpose(score)
