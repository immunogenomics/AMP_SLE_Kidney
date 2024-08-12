import scanpy as sc
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import scipy.stats as st
import cna
import multianndata as mad
np.random.seed(0)


def run_cna(variable, cna_obj, covs = None):
    if covs is not None:
        cna_res = cna.tl._association.association(cna_obj, #dataset 
                                                  cna_obj.samplem[variable], #phenotype
            #                                       batches=d.samplem.processing_batch, #batches
                                                  covs = cna_obj.samplem[covs],
                                                  Nnull=10000, # number of null permutations to do (defaults to only 1e3)
                                                  ks=[1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20] # I asked the method to consider up to 10 PCs because
                                                                            #it chose the max number of PCs it considered the default set of [1,2,3,4]
                                                 )
    else:
        cna_res = cna.tl._association.association(cna_obj, #dataset 
                                                  cna_obj.samplem[variable], #phenotype
            #                                       batches=d.samplem.processing_batch, #batches
                                                 # covs = d.samplem[covs]
                                                  Nnull=10000, # number of null permutations to do (defaults to only 1e3)
                                                  ks=[1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20] # I asked the method to consider up to 10 PCs because
                                                                            #it chose the max number of PCs it considered the default set of [1,2,3,4]
                                                 )
        
    return cna_res
def cna_test(meta_df, harmony_df, umap_df, variable, covars = None, ret_cna_obj = False):
    print(meta_df.shape)
    print(harmony_df.shape)
    print(umap_df.shape)
    cna_object = mad.MultiAnnData(X=harmony_df, obs=meta_df, sampleid="sample")
    if covars is not None:   
        cna_object.obs_to_sample(covars + [variable])
        umap_df.index = cna_object.obs.index
        cna_object.obsm['X_umap'] = umap_df
        cna.pp.knn(cna_object)
        np.random.seed(0) 
        final_res = run_cna(variable, cna_object, covs = covars)
    else:
        cna_object.obs_to_sample([variable])
        umap_df.index = cna_object.obs.index
        cna_object.obsm['X_umap'] = umap_df
        cna.pp.knn(cna_object)
        np.random.seed(0) 
        final_res = run_cna(variable, cna_object, covs = None)
    print(final_res.k)
    print(final_res.ks)
    print('p =', final_res.p, ',', final_res.k, 'PCs used')
    print('total r^2 between top {} NAM PCs and outcome is {:.2f}'.format(final_res.k, final_res.r2))
    if ret_cna_obj:
        return cna_object, final_res
    else:
        return final_res