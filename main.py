import scanpy as sc
from GeneClust import select

test_adata = sc.datasets.pbmc3k()
sc.pp.filter_cells(test_adata, min_genes=10)
sc.pp.filter_genes(test_adata, min_cells=10)
test_adata.X = test_adata.X.toarray()
test_adata.raw = test_adata
sc.pp.normalize_total(test_adata)
sc.pp.log1p(test_adata)
test_adata.uns['data_name'] = 'pbmc3k'

selected = select(test_adata, n_selected_genes=500, dr_method='pca', n_comps=50, similarity='spearman',
                  clustering='gmm', n_clusters=200, in_cluster_score='m3drop', return_genes=True)
print(selected)

# TODO: umap

