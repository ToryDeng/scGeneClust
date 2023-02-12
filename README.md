![Github license](https://img.shields.io/github/license/ToryDeng/scGeneClust)
![Github language](https://img.shields.io/pypi/pyversions/GeneClust)
![Github version](https://img.shields.io/pypi/v/GeneClust)

# **GeneClust**: a cofunctional grouping-based approach for non-redundant feature gene selection in unannotated single-cell RNA-seq

GeneClust is a computational feature selection method for scRNA-seq cell clustering. GeneClust groups genes into clusters from which genes are evaluated and selected with the aim of maximizing *relevance*, minimizing *redundancy* and preserving *complementarity*.
![image](https://github.com/ToryDeng/scGeneClust/blob/main/docs/images/workflow.png?raw=true)

## Dependencies

- anndata>=0.8.0
- numpy>=1.21.6
- setuptools>=59.5.0
- scanpy>=1.9.1
- scipy>=1.9.3
- loguru>=0.6.0
- hdbscan>=0.8.29
- sklearn>=0.0.post2
- scikit-learn>=1.2.0
- igraph>=0.10.2
- leidenalg>=0.9.1
- pandas>=1.5.2
- SpaGCN>=1.2.5
- squidpy>=1.2.2
- torch>=1.13.1
- opencv-python>=4.6.0

## Installation

1. **PyPI**

You can directly install the package from PyPI.

```
pip3 install GeneClust
```

2. **Github**

Also, You can download the package from GitHub and install it locally:

```
git clone https://github.com/ToryDeng/scGeneClust.git
cd scGeneClust/
python3 setup.py install --user
```

## Two Versions of GeneClust

|  **Version**  | **Usage Scenarios**                                                                                       |
| :--------------: | ----------------------------------------------------------------------------------------------------------- |
|  GeneClust-ps  | 1. The number of cells is small (e.g., several thousand) 2. Cell clustering performance is more important |
| GeneClust-fast | 1. The number of cells is large (e.g., over 50,000) 2. Computational efficiency is more important         |

## Tutorial

For the step-by-step tutorial, please refer to the notebook:

https://github.com/ToryDeng/scGeneClust/blob/main/notebooks/tutorial_scRNA-seq.ipynb


## Reproducibility
To reproduce the results and figures presented in our paper, please go to https://github.com/ToryDeng/scGeneClust/tree/main/figures.

## Citation
```bibtex
@article{10.1093/bib/bbad042,
    author = {Deng, Tao and Chen, Siyu and Zhang, Ying and Xu, Yuanbin and Feng, Da and Wu, Hao and Sun, Xiaobo},
    title = "{A cofunctional grouping-based approach for non-redundant feature gene selection in unannotated single-cell RNA-seq analysis}",
    journal = {Briefings in Bioinformatics},
    year = {2023},
    month = {02},
    issn = {1477-4054},
    doi = {10.1093/bib/bbad042},
    url = {https://doi.org/10.1093/bib/bbad042},
    note = {bbad042},
    eprint = {https://academic.oup.com/bib/advance-article-pdf/doi/10.1093/bib/bbad042/49126753/bbad042.pdf},
}
```