# -*- coding: utf-8 -*-
# @Time : 2022/6/6 17:19
# @Author : Tory Deng
# @File : setup.py
# @Software: PyCharm
from setuptools import find_packages
from setuptools import setup

with open('README.md') as f:
    LONG_DESCRIPTION = f.read()

setup(
    name='GeneClust',
    version='1.0.0',
    description='A cofunctional grouping-based approach for non-redundant feature gene selection in '
                'unannotated single-cell RNA-seq',
    long_description=LONG_DESCRIPTION,
    long_description_content_type='text/markdown',
    author='Tao Deng',
    author_email='taodeng@link.cuhk.edu.cn',
    url='https://github.com/ToryDeng/scGeneClust',
    license='GPL v3',
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: GNU General Public License v3 (GPLv3)",
        "Operating System :: OS Independent"
    ],
    python_requires=">=3.9",
    install_requires=[
        'anndata>=0.8.0',
        'numpy>=1.21.6',
        'setuptools>=59.5.0',
        'scanpy>=1.9.1',
        'scipy>=1.9.3',
        'loguru>=0.6.0',
        'hdbscan>=0.8.29',
        'sklearn>=0.0.post2',
        'scikit-learn>=1.2.0',
        'igraph>=0.10.2',
        'leidenalg>=0.9.1',
        'pandas>=1.5.2',
        'SpaGCN>=1.2.5',
        'squidpy>=1.2.2',
        'torch>=1.13.1',
        'opencv-python>=4.6.0',
    ],
    packages=find_packages(exclude=('tests', 'figures', 'data', 'docs', 'notebooks')),
    zip_safe=False,
)
