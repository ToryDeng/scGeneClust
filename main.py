import numpy as np
import scanpy as sc
import numba
from GeneClust import select, load_example_adata
from sigclust import SigClust


test_adata = load_example_adata()
# selected = select(test_adata, n_selected_genes=500, dr_method='pca', n_comps=50, distance='rho_p',
#                   clustering='ms', in_cluster_score='center', inter_cluster_score='silhouette',
#                   return_genes=True)

# print(test_adata)
