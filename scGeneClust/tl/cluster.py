
import numpy as np
from sklearn.feature_selection import mutual_info_regression
import networkx as nx
import anndata as ad
from typing import Optional, Literal
from joblib import Parallel, delayed


def clustering(
        adata: ad.AnnData,
        scale: int=2000,
        random_stat: int=42):
    
    rlv_mtx = adata.layers['X_gene_log'][adata.obs['highly_confident'] == True][:,adata.var['relevant']]  

    #Calculate mi between genes to get an upper triangular matrix
    mi_mtx = np.zeros(shape=(rlv_mtx.shape[1],rlv_mtx.shape[1]))
    res = Parallel(n_jobs=-1)(delayed(mutual_info_regression)(rlv_mtx[:,i+1:], rlv_mtx[:,i]) for i in range(rlv_mtx.shape[1]-1))
    for i in range(len(res)):
        mi_mtx[i, i+1:] = res[i]

    #Constructe MST 
    G = nx.Graph(mi_mtx) 
    MST = nx.minimum_spanning_tree(G,weight="weight", algorithm="prim")  

    #Prune MST: mi < max{comp(mode_1,node_2),min{rlv(node_1),rlv(node_2)}}
    results = Parallel(n_jobs=-1)(delayed(prune)(MST, adata, scale) for i in range(1))
    MST_pruned=results[0][0]   

    #Cluster:find subtrees
    clusters = list(nx.connected_components(MST_pruned)) 
    adata.var['cluster'] = None
    data_reset = adata.var.reset_index(inplace=False)
    for i in range(len(clusters)):
        if len(clusters[i]) == 1:
            loc = data_reset[data_reset['index'] == adata.var[adata.var['relevant'] == True].iloc[list(clusters[i])[0]].name].index[0]
            adata.var['cluster'][loc] = -1
        else:
            for j in clusters[i]:
                loc = data_reset[data_reset['index'] == adata.var[adata.var['relevant'] == True].iloc[j].name].index[0]
                adata.var['cluster'][loc] = i

    #Find representative gene in each non single cluster and single clusters
        #Find the mininum score of representative gene in all non single clusters
    grouped = adata.var[(adata.var['relevant'] == True)&(adata.var['cluster'] != -1)].groupby(by='cluster')['relevance']
    thres = min(grouped.nlargest(1))
        #Single gene index
    single_repre__gene_index = adata.var[(adata.var['relevant'] == True)&(adata.var['cluster'] == -1)&(adata.var['relevance']>thres)].index
        #Non single gene index
    nonsingle_rep_gene_index = [grouped.nlargest(1).index[i][1] for i in range(len(grouped.nlargest(1)))]

    selected_index = []
    selected_index.extend(single_repre__gene_index)
    selected_index.extend(nonsingle_rep_gene_index)
    # adata.var['rep'] = False
    # adata.var['rep'][single_repre__gene_index] = True
    # adata.var['rep'][nonsingle_rep_gene_index] = True

    selected_genes = selected_index
    return selected_genes

    
def complementarity(
                adata: ad.AnnData,
                node1: int, 
                node2: int, 
                random_stat: int):

    data_reset = adata.var.reset_index()
    cmi = 0
    for clus in np.unique(adata.obs['cluster']):
        if clus != -1:
            node1_loc = data_reset[data_reset['index'] == adata.var[adata.var['relevant'] == True].iloc[node1].name].index[0]
            node2_loc = data_reset[data_reset['index'] == adata.var[adata.var['relevant'] == True].iloc[node2].name].index[0]
            # rlv = Parallel(n_jobs=-1)(delayed(mutual_info_regression)(adata.layers['X_gene_log'][adata.obs['cluster'] == clus][:,node1_loc].reshape(-1,1),adata.layers['X_gene_log'][adata.obs['cluster'] == clus][:,node2_loc]) for i in range(1))
            # cmi += rlv[0][0]*(len(adata.obs[adata.obs['cluster'] == clus])/len(adata.obs[adata.obs['highly_confident'] == True]))
            rlv = mutual_info_regression(adata.layers['X_gene_log'][adata.obs['cluster'] == clus][:,node1_loc].reshape(-1,1),adata.layers['X_gene_log'][adata.obs['cluster'] == clus][:,node2_loc],random_state=random_stat)
            cmi += rlv*(len(adata.obs[adata.obs['cluster'] == clus])/len(adata.obs[adata.obs['highly_confident'] == True]))
    return cmi

def prune(MST, adata, scale):
    MST_pruned = MST.copy()
    complms = []
    class_rlvs = []
    weights = []
    for node_1, node_2, w in MST.edges(data=True):
        class_rlv = min(adata.var[adata.var['relevant'] == True].iloc[node_1]['relevance'], adata.var[adata.var['relevant'] == True].iloc[node_2]['relevance'])
        complm = complementarity(adata, node_1, node_2, 42) 
        class_rlvs.append(class_rlv)
        complms.append(complm[0])
        weights.append(w['weight'])
        if w['weight']*(scale) < max(class_rlv, complm):
            MST_pruned.remove_edge(node_1, node_2)
    return MST_pruned, complms, class_rlvs, weights