import os
import pickle
import scanpy as sc
from scipy import sparse
import numpy as np
from sklearn.model_selection import train_test_split
from anndata import AnnData

def pklthat(gdata, fname, fpath): 
    with open(os.path.join(fpath,fname),'wb') as f :
        pickle.dump(gdata, f, protocol=pickle.HIGHEST_PROTOCOL)
        f.close()

def graph_pp(AnnData:AnnData, num_comps:int,num_neighbors:int,bbknn=False)->AnnData:
    """
    Calculate neighbors graph
    
    param Anndata: anndata object
    param num_comps: numper of principal components
    param num_neighbors: number of neighbors
    
    return Anndata: anndata object
    """
    sc.tl.pca(AnnData, n_comps=num_comps)
    if bbknn:
        sc.external.pp.bbknn(AnnData) # use default params
    else:
        sc.pp.neighbors(AnnData, n_pcs=num_comps, n_neighbors=num_neighbors)
    return AnnData

def dictthat(AnnData:AnnData, gene_ranger=True)->dict:
    """
    Create object with graph representation of data
    
    param AnnData: Anndata object
    param gene_ranger: bool 
    
    return gdata: dict which contain graph representation of Anndata object
    """
    if gene_ranger:
        minimum = AnnData.X.min(axis=0)
        maximum = AnnData.X.max(axis=0)
        if sparse.issparse(AnnData.X):
            num = AnnData.X - minimum.todense()
            denom =  (maximum - minimum).todense()
        else:
            num = AnnData.X - minimum
            denom =  maximum - minimum
        xhat = np.divide(num, denom, out=np.zeros_like(num), where=denom!=0) 
    else:
        xhat = (AnnData.X - AnnData.X.min()) / (AnnData.X.max() - AnnData.X.min())
    gdata = {'X':xhat,
             'adj':AnnData.uns['neighbors']['connectivities']+sparse.diags([1]*AnnData.shape[0], format='csr'),
             'feature_names':AnnData.var_names.to_list()}
    gdata['cell_id'] = AnnData.obs.index.to_list()
    for col in AnnData.obs.columns:
        gdata[col] = AnnData.obs[col].to_list()
    
    return gdata


def patients_genotype_check(adata_obs,gene_col:str,patient_col:str):
    for gene in set(adata_obs[gene_col]):
        print('gene: ', gene, 'num patients: ',len(set(adata_obs[adata_obs[gene_col]==gene][patient_col])))
        print('genotype samples',np.unique(adata_obs[adata_obs[gene_col]==gene][gene_col],return_counts = True))