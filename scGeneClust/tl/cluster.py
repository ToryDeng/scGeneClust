import numpy as np
from sklearn.feature_selection import mutual_info_regression
import networkx as nx
import anndata as ad
from typing import Optional, Literal
from joblib import Parallel, delayed
from scGeneClust.tl import handle_single_gene_cluster
from loguru import logger


def clustering(
        adata: ad.AnnData,
        scale: int = 2000,
        random_stat: int = 42
):
    logger.info(f"Start to cluster genes...")
    # Calculate mi between genes to get an upper triangular matrix
    mi_mtx = np.zeros(shape=(adata.layers['X_gene_log'].shape[1], adata.layers['X_gene_log'].shape[1]))
    res = Parallel(n_jobs=80)(delayed(mutual_info_regression)
                              (adata.layers['X_gene_log'][:, i + 1:], adata.layers['X_gene_log'][:, i])
                              for i in range(adata.layers['X_gene_log'].shape[1] - 1))
    for i in range(len(res)):
        mi_mtx[i, i + 1:] = res[i]

    # Construct MST
    G = nx.Graph(mi_mtx)
    MST = nx.minimum_spanning_tree(G, weight="weight", algorithm="prim")
    logger.debug("MST constructed!")

    # Prune MST: mi < max{comp(mode_1,node_2),min{rlv(node_1),rlv(node_2)}}
    node_pairs = Parallel(n_jobs=80)(delayed(prune)(adata, node1, node2, w['weight'], scale, random_stat)
                                     for node1, node2, w in MST.edges(data=True))
    MST_pruned = MST.copy()
    for node_pair in node_pairs:
        if node_pair is not None:
            MST_pruned.remove_edge(node_pair[0], node_pair[1])
    logger.debug("MST pruned!")

    # Cluster:find subtrees
    clusters = list(nx.connected_components(MST_pruned))
    logger.debug(f"number of clusters: {len(clusters)}")
    adata.var['cluster'] = None
    for i in range(len(clusters)):
        if len(clusters[i]) == 1:
            adata.var['cluster'][list(clusters[i])[0]] = -1
        else:
            for j in clusters[i]:
                adata.var['cluster'][j] = i

    # Find representative genes
    adata.var['representative'] = False
    grouped = adata.var[adata.var['cluster'] != -1].groupby(by='cluster')['score']
    for i in range(len(grouped.nlargest(1))):
        adata.var.loc[grouped.nlargest(1).index[i][1], 'representative'] = True

    logger.debug(f"number of rep genes in non single cluster: {len(adata.var[adata.var['representative']])}")
    if len(adata.var[adata.var['cluster'] == -1]) > 0:
        handle_single_gene_cluster(adata, mode='hc', random_stat=random_stat)

    logger.info(f"Gene clustering finished...")


def prune(adata, node1, node2, weight, scale, random_stat):
    class_rlv = min(adata.var['score'][node1], adata.var['score'][node2])
    complm = complementarity(adata, node1, node2, random_stat)
    if weight * scale < max(class_rlv, complm):
        return node1, node2
    else:
        return None


def complementarity(
        adata: ad.AnnData,
        node1: int,
        node2: int,
        random_stat: int
):
    cmi = 0
    for clus in np.unique(adata.obs['cluster']):
        rlv = mutual_info_regression(adata.layers['X_gene_log'][adata.obs['cluster'] == clus][:, node1].reshape(-1, 1),
                                     adata.layers['X_gene_log'][adata.obs['cluster'] == clus][:, node2],
                                     random_state=random_stat)
        cmi += rlv * (len(adata.obs[adata.obs['cluster'] == clus]) / len(adata.obs))
    return cmi
