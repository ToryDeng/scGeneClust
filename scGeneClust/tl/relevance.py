import anndata as ad
from typing import Optional
from sklearn.feature_selection import mutual_info_classif
from loguru import logger

def find_relevant_gene(
        adata: ad.AnnData,
        rlv_threshold: float = 0.01,
        random_stat: int = 42):

    #find relevant gene according to mutual information with cluster label based on high confident cells
    logger.info(f"Start to find relevant genes...")
    relevance = mutual_info_classif(adata.layers['X_gene_log'][adata.obs['highly_confident'] == True],adata.obs.cluster[adata.obs['highly_confident'] == True],discrete_features=True,random_state=random_stat)
    adata.var['relevant'] = relevance > rlv_threshold
    adata.var['relevance'] = relevance
    logger.info(f"Relevant gene detection finished!")




