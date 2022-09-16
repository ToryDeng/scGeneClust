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
    version='0.0.1',
    description='Cofunctional grouping-based feature gene selection for unsupervised scRNA-seq clustering',
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
    python_requires=">=3.8",
    install_requires=[
        'pandas>=1.4.2',
        'anndata>=0.8.0',
        'numpy>=1.21.5',
        'loguru>=0.6.0',
        'sklearn>=0.0',
        'scikit-learn>=1.1.1',
        'scanpy>=1.9.1',
        'scipy>=1.7.3',
        'leidenalg>=0.8.9'
    ],
    packages=find_packages(exclude=('tests', 'figures')),
    zip_safe=False,
)
