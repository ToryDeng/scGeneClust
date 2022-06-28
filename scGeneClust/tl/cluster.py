import numpy as np
from sklearn.feature_selection import mutual_info_regression
import networkx as nx
import anndata as ad
from itertools import combinations
from functools import partial
from joblib import Parallel, delayed
from typing import Optional
from multiprocessing.pool import Pool
from multiprocessing import cpu_count
from scipy.spatial.distance import squareform
from scGeneClust.tl import handle_single_gene_cluster
from loguru import logger


def clustering(
        adata: ad.AnnData,
        scale: int = 2000,
        random_stat: int = 42
):
    logger.info(f"Start to cluster genes...")

    # calculate mi between genes to get an upper triangular matrix
    pool = Pool(processes=cpu_count() - 1)
    partial_cal_mi = partial(cal_mi, data=adata.layers['X_gene_log'], random_stat=random_stat)
    mi_mtx = squareform(pool.starmap(partial_cal_mi, combinations(range(adata.n_vars), 2)))

    # construct MST
    G = nx.Graph(mi_mtx)
    MST = nx.minimum_spanning_tree(G, weight="weight", algorithm="prim")
    logger.debug("MST constructed!")

    # prune MST: mi < max{comp(mode_1,node_2),min{rlv(node_1),rlv(node_2)}}
    node_pairs = pool.starmap(partial(prune, adata=adata, scale=scale, random_stat=random_stat), MST.edges(data=True))
    for node_pair in node_pairs:
        if node_pair is not None:
            MST.remove_edge(node_pair[0], node_pair[1])
    logger.debug("MST pruned!")

    # cluster genes by finding subtrees
    clusters = list(nx.connected_components(MST))
    adata.var['cluster'] = np.hstack([np.full((len(cluster),), fill_value=i) for i, cluster in enumerate(clusters)])
    logger.info(f"Gene clustering finished...")


def cal_mi(i: int, j: int, data: np.ndarray, random_stat: Optional[int]):
    return mutual_info_regression(
        data[:, i].reshape(-1, 1), data[:, j], discrete_features=False, random_state=random_stat
    )[0]


def prune(node1: int, node2: int, w, adata, scale, random_stat):
    class_rlv = min(adata.var['score'][node1], adata.var['score'][node2])
    complm = complementarity(adata, node1, node2, random_stat)
    return (node1, node2) if w['weight'] * scale < max(class_rlv, complm) else None


def complementarity(adata: ad.AnnData, n1: int, n2: int, random_stat: Optional[int]):
    cmi = 0
    data = adata.layers['X_gene_log']
    for clus in np.unique(adata.obs['cluster']):
        clus_mask = adata.obs['cluster'] == clus
        rlv = mutual_info_regression(
            data[clus_mask, n1].reshape(-1, 1), data[clus_mask, n2], discrete_features=False, random_state=random_stat
        )
        cmi += rlv * clus_mask.sum() / adata.n_obs
    return cmi
