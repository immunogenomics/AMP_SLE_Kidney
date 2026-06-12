import multianndata as mad
import anndata as ad
import pandas as pd
import cna
import sys
import scanpy as sc
import matplotlib.pyplot as plt
import random
import numpy as np

np.random.seed(0) 
random.seed(0)

name=sys.argv[1]
sampleType=sys.argv[2]

if sampleType == "tissue":
    folderName = "cna_new"
    harmony = pd.read_csv(('/data/srlab/ssg34/SLE_kidney_v2/data/' +  folderName + '/' + name + '/activity/sc_harmony.csv'))
    meta = pd.read_csv(('/data/srlab/ssg34/SLE_kidney_v2/data/' +  folderName + '/' + name + '/activity/sc_meta.csv'))
    umap = pd.read_csv(('/data/srlab/ssg34/SLE_kidney_v2/data/' +  folderName + '/' + name + '/activity/sc_umap.csv'))

    d = mad.MultiAnnData(X=harmony, obs=meta, sampleid="sample")
else:
    folderName = "pbmc"
    harmony = pd.read_csv(('/data/srlab/ssg34/SLE_kidney_v2/data/' +  folderName + '/' + name + '/activity/harmony.csv'))
    meta = pd.read_csv(('/data/srlab/ssg34/SLE_kidney_v2/data/' +  folderName + '/' + name + '/activity/meta.csv'))
    umap = pd.read_csv(('/data/srlab/ssg34/SLE_kidney_v2/data/' +  folderName + '/' + name + '/activity/umap.csv'))
    
    d = mad.MultiAnnData(X=harmony, obs=meta, sampleid="Unified_Visit")

d.obs_to_sample(['Final_Chronicity', 'Final_Activity', 'First_biop', 'Responder_Status', 
       'Race_[W]', 'Final_ISN_[III]', 'Final_ISN_[III][V]', 'Final_ISN_[IV]',
       'Final_ISN_[IV][V]', 'Final_ISN_[V]'])

umap.index = d.obs.index
d.obsm['X_umap'] = umap

cna.pp.knn(d)

def run_cna(variable, cna_obj, fileName, covs = None):
    if covs is not None:
        cna_res = cna.tl._association.association(cna_obj, #dataset 
                                                  cna_obj.samplem[variable], #phenotype
                                                  covs = cna_obj.samplem[covs],
                                                  Nnull=10000, # number of null permutations to do (defaults to only 1e3)
                                                  ks=[1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20] # I asked the method to consider up to 10 PCs because
                                                                            #it chose the max number of PCs it considered the default set of [1,2,3,4]
                                                 )
    else:
        cna_res = cna.tl._association.association(cna_obj, #dataset 
                                                  cna_obj.samplem[variable], #phenotype
                                                  Nnull=10000, # number of null permutations to do (defaults to only 1e3)
                                                  ks=[1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20] # I asked the method to consider up to 10 PCs because
                                                                            #it chose the max number of PCs it considered the default set of [1,2,3,4]
                                                 )
    print("p = " + str(cna_res.p) + " with a minimum FDR of " + str(min(cna_res.fdrs.fdr)))
    np.savetxt("cna_results/" + name + "_" + sampleType + "/" + variable + "_" + fileName + "_ncorr.csv", 
                cna_res.ncorrs, delimiter=",")
    np.savetxt("cna_results/" + name + "_" + sampleType + "/" + variable + "_" + fileName + "_fdrs.csv", 
                cna_res.fdrs, delimiter=",")
    return cna_res

print("Sanity check!")
activityCNA = run_cna('Final_Chronicity', d, "justChecking")

print("Testing Activity!")

activityCNA_responder = run_cna('Final_Activity', d, "noISN", ['First_biop', 'Final_Chronicity', 'Responder_Status'] )
activityCNA_responder = run_cna('Final_Activity', d, "withISN_IV_V", ['First_biop', 'Final_Chronicity', 'Responder_Status', "Final_ISN_[IV][V]"] )
activityCNA_responder = run_cna('Final_Activity', d, "withISN_V", ['First_biop', 'Final_Chronicity', 'Responder_Status', "Final_ISN_[V]"] )
