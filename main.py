import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from sklearn.metrics.cluster import contingency_matrix

from pagest import pagest
from pagest.utils import load_example_adata

adata = load_example_adata(min_cells=2)
genes = pagest(adata, n_features=500, mode='one-way', verbose=2)
adata.write(f"cache/{adata.uns['data_name']}_pagest.h5ad")
contingency_df = pd.DataFrame(contingency_matrix(adata.obs.celltype, adata.obs.cluster),
                              index=np.unique(adata.obs.celltype), columns=np.unique(adata.obs.cluster))
plt.figure(figsize=(15, 9))
sns.heatmap(contingency_df, cmap='YlGnBu', annot=True, square=False, fmt='d')
plt.savefig(f"figures/{adata.uns['data_name']}_heatmap.svg", dpi=300)

