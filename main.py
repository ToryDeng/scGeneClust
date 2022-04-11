import scanpy as sc
from GeneClust import select


# 1000 cells, 5000 genes after log-norm

test_adata = sc.datasets.pbmc3k()
test_adata.X = test_adata.X.toarray()
test_adata.raw = test_adata
sc.pp.recipe_seurat(test_adata)

selected = select(test_adata, n_selected_genes=201, dr_method='glm-pca', n_comps=50, similarity='pearson', n_clusters=200, return_genes=True)
print(selected)
