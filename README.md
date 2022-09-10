![Github license](https://img.shields.io/github/license/ToryDeng/scGeneClust)
![Github language](https://img.shields.io/pypi/pyversions/GeneClust)
![Github version](https://img.shields.io/pypi/v/GeneClust)

# **GeneClust**: cofunctional grouping-based feature gene selection for unsupervised scRNA-seq clustering
GeneClust is a computational feature selection method for scRNA-seq cell clustering. GeneClust groups genes into clusters from which genes are evaluated and selected with the aim of maximizing *relevance*, minimizing *redundancy* and preserving *complementarity*. 
![image](https://github.com/ToryDeng/scGeneClust/blob/main/docs/images/workflow.png)
## Dependencies
- numpy>=1.21.5
- pandas>=1.4.2
- anndata>=0.8.0
- setuptools>=59.5.0
- loguru>=0.6.0
- sklearn>=0.0
- scikit-learn>=1.1.1
- scanpy>=1.9.1
- scipy>=1.7.3
- leidenalg>=0.8.9
## Installation
1. **PyPI**

You can directly install the package from PyPI.
```
pip3 install GeneClust
```

2. **Github**

Also, You can download the package from Github and install it locally:
```
git clone https://github.com/ToryDeng/scGeneClust.git
cd scGeneClust/
python3 setup.py install --user
```
## Two Versions of GeneClust
| **Version** | **Usage Scenarios** |
|   :----:   |   --------   |
|  GeneClust-ps | 1. Number of cells is small (e.g., several thousand) <br> 2. Cell clustering performance is more important  |
|  GeneClust-fast   |    1. Number of cells is large (e.g., over 50,000) <br> 2. Computational efficiency is more important   |
## Example Code
```python
from scGeneClust.utils import load_PBMC3k
from scGeneClust import scGeneClust

# load the PBMC3k dataset
raw_adata = load_PBMC3k()
# GeneClust-fast
selected_genes = scGeneClust(raw_adata, version='fast', n_gene_clusters=200, random_stat=2022, verbosity=2)
# GeneClust-ps
selected_genes = scGeneClust(raw_adata, version='ps', n_cell_clusters=7, scale=1000, top_percent_relevance=5, random_stat=2022, verbosity=2)
```
GeneClust expects raw counts. The output is an ndarray of selected features, which can be used in the downstream cell clustering analysis.



