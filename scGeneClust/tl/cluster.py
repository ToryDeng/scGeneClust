from functools import partial
from itertools import combinations
from multiprocessing import cpu_count
from multiprocessing.pool import Pool
from typing import Optional

import anndata as ad
import networkx as nx
import numpy as np
import scGeneClust.tl as tl
from loguru import logger
from scipy.spatial.distance import squareform
from sklearn.cluster import MiniBatchKMeans
from sklearn.feature_selection import mutual_info_regression


def gene_clustering_mbkmeans(
        adata: ad.AnnData,
        n_gene_clusters: int,
        random_stat: Optional[int]
):
    km = MiniBatchKMeans(n_clusters=n_gene_clusters, batch_size=1024, random_state=random_stat)
    adata.var['cluster'] = km.fit_predict(adata.varm['pca'])  # gene gene_clustering_graph
    adata.var['score'] = tl.compute_gene_closeness(adata, km.cluster_centers_)


def gene_clustering_graph(
        adata: ad.AnnData,
        scale: int,
        random_stat: Optional[int]
):
    logger.info(f"Start to cluster genes...")

    # calculate mi between genes to get an upper triangular matrix
    pool = Pool(processes=cpu_count() - 1)  #
    partial_cal_mi = partial(cal_mi, data=adata.layers['X_gene_log'], random_stat=random_stat)
    adata.varp['redundancy'] = squareform(pool.starmap(partial_cal_mi, combinations(range(adata.n_vars), 2)))

    # construct MST
    G = nx.Graph(adata.varp['redundancy'])
    MST = nx.minimum_spanning_tree(G, weight="weight", algorithm="prim")
    logger.debug("MST constructed!")

    # prune MST: mi < max{comp(mode_1,node_2),min{rlv(node_1),rlv(node_2)}}
    partial_prune = partial(prune,
                            data=adata.layers['X_gene_log'], clusters=adata.obs['cluster'], rlv=adata.var['score'],
                            scale=scale, random_stat=random_stat)
    node_pairs = pool.starmap(partial_prune, MST.edges(data=True))
    for node_pair in node_pairs:
        if node_pair is not None:
            MST.remove_edge(node_pair[0], node_pair[1])
    logger.debug("MST pruned!")

    # cluster genes by finding subtrees
    clusters = list(nx.connected_components(MST))
    adata.var['cluster'] = np.hstack([np.full((len(cluster),), fill_value=i) for i, cluster in enumerate(clusters)])
    logger.info(f"Gene gene_clustering_graph finished...")


def cal_mi(i: int, j: int, data: np.ndarray, random_stat: Optional[int]):
    return mutual_info_regression(
        data[:, i].reshape(-1, 1), data[:, j], discrete_features=False, random_state=random_stat
    )[0]


def prune(node1: int, node2: int, w, data, clusters, rlv, scale, random_stat):
    class_rlv = min(rlv[node1], rlv[node2])
    complm = cal_complementarity(data, clusters, node1, node2, random_stat)
    return (node1, node2) if w['weight'] * scale < max(class_rlv, complm) else None


def cal_complementarity(data, clusters, n1: int, n2: int, random_stat: Optional[int]):
    cmi = 0
    for clus in np.unique(clusters):
        clus_mask = clusters == clus
        rlv = mutual_info_regression(
            data[clus_mask, n1].reshape(-1, 1), data[clus_mask, n2], discrete_features=False, random_state=random_stat
        )
        cmi += rlv * clus_mask.sum() / data.shape[0]
    return cmi
