'''
Description: Similarity measure based on Bayes Correlation
Author: Ying Zhang
Reference: Bayesian Correlation is a robust gene similarity measure for single cell RNA-seq data
'''


import pandas as pd
import numpy as np
def bayes_corr(data: pd.DataFrame) -> pd.DataFrame:

    """
    similarity measure using Bayesian correlation 
    :param data: columns: genes 
    :return: The similarity matrix, in which each entry is the similarity between two genes
    """

    data = data.T
    nrowsX = data.shape[0]
    ncolsX = data.shape[1]
    
    alpha0 = [1/nrowsX]*ncolsX
    beta0 = [1-x for x in alpha0]
    cs = data.sum(axis = 0).tolist()
    alphas  = np.array(data) + alpha0
    betas = np.array(beta0*nrowsX).reshape(nrowsX,-1)+np.array(cs*nrowsX).reshape(nrowsX,-1)-np.array(data)
    alphasPLUSbetas  = alphas + betas
    alp_alpPLUSbeta = alphas/alphasPLUSbetas 
    Psi = alp_alpPLUSbeta - np.array([x/ncolsX for x in alp_alpPLUSbeta.sum(axis = 1)]*ncolsX).reshape(-1,nrowsX).T
    var_vec = ((((alphas*betas)/((alphasPLUSbetas**2)*(alphasPLUSbetas+1))).sum(axis = 1)+(Psi**2).sum(axis = 1))/ncolsX).reshape(nrowsX,1)
    cov_mtrx = np.dot(Psi,Psi.T)/ncolsX
    Bcorrvals = cov_mtrx/np.sqrt(np.dot(var_vec,var_vec.T))
    Bcorrvals[np.diag_indices_from(Bcorrvals)] = 1
    return pd.DataFrame(Bcorrvals)




    






