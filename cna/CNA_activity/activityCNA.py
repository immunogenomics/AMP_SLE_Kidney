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

harmony = pd.read_csv(('/data/srlab/ssg34/SLE_kidney_v2/data/cna_new/' + name + '/activity/sc_harmony.csv'))
meta = pd.read_csv(('/data/srlab/ssg34/SLE_kidney_v2/data/cna_new/' + name + '/activity/sc_meta.csv'))
umap = pd.read_csv(('/data/srlab/ssg34/SLE_kidney_v2/data/cna_new/' + name + '/activity/sc_umap.csv'))

d = mad.MultiAnnData(X=harmony, obs=meta, sampleid="sample")

d.obs_to_sample(["Final_Activity", "Final_Chronicity", "First_biop", "Responder_Status", 
                 "Final_ISN_[III]", "Final_ISN_[III][V]", "Final_ISN_[IV]", "Final_ISN_[IV][V]", "Final_ISN_[V]",
                 "Final_Site_Cincinnati", "Final_Site_Einstein", "Final_Site_JHU", 
                 "Final_Site_MUSC", "Final_Site_Northwell", "Final_Site_NYU", "Final_Site_Rochester", 
                 "Final_Site_UCLA", "Final_Site_UCSF", "Pred_use"])

umap.index = d.obs.index
d.obsm['X_umap'] = umap

cna.pp.knn(d)

with open((name +  "_pValues.csv"), 'a') as file:
    # Append the line of text
    file.write("name,pValue\n")

def run_cna(variable, cna_obj, fileName, covs = None):
    if covs is not None:
        cna_res = cna.tl._association.association(cna_obj, #dataset 
                                                  cna_obj.samplem[variable], #phenotype
                                                  covs = cna_obj.samplem[covs],
                                                  Nnull=100000, # number of null permutations to do (defaults to only 1e3)
                                                  ks=[1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20] # I asked the method to consider up to 10 PCs because
                                                                            #it chose the max number of PCs it considered the default set of [1,2,3,4]
                                                 )
    else:
        cna_res = cna.tl._association.association(cna_obj, #dataset 
                                                  cna_obj.samplem[variable], #phenotype
                                                  Nnull=100000, # number of null permutations to do (defaults to only 1e3)
                                                  ks=[1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20] # I asked the method to consider up to 10 PCs because
                                                                            #it chose the max number of PCs it considered the default set of [1,2,3,4]
                                                 )
    with open((name +  "_pValues.csv"), 'a') as file:
        # Append the line of text
        file.write(fileName + "," + str(cna_res.p) + '\n')
    np.savetxt("cna_results/" + name + "/" + variable + "_" + fileName + "_ncorr.csv", 
                cna_res.ncorrs, delimiter=",")
    np.savetxt("cna_results/" + name + "/" + variable + "_" + fileName + "_fdrs.csv", 
                cna_res.fdrs, delimiter=",")
    return cna_res

chronicityCNA = run_cna('Final_Activity', d, "None")

chronicityCNA_SiteFirstBiop = run_cna('Final_Activity', d, "SiteFirstBiopsy", 
                                              ["Final_Site_Einstein", "Final_Site_JHU", 
                                              "Final_Site_MUSC", "Final_Site_NYU", "Final_Site_Rochester", 
                                              "Final_Site_UCSF", 'First_biop'])

chronicityCNA_SiteFirstBiopChron = run_cna('Final_Activity', d, "SiteFBChron", 
                                              ["Final_Site_Einstein", "Final_Site_JHU", 
                                              "Final_Site_MUSC", "Final_Site_NYU", "Final_Site_Rochester", 
                                              "Final_Site_UCSF", 'First_biop', 'Final_Chronicity'])

chronicityCNA_ISN_V = run_cna('Final_Activity', d, "SiteFBChron_ISN_V", 
                                              ["Final_Site_Einstein", "Final_Site_JHU", 
                                              "Final_Site_MUSC", "Final_Site_NYU", "Final_Site_Rochester", 
                                              "Final_Site_UCSF", 'First_biop', 'Final_Chronicity', 
                                              "Final_ISN_[V]"])

chronicityCNA_ISN_IV = run_cna('Final_Activity', d, "SiteFBChron_ISN_IV", 
                                              ["Final_Site_Einstein", "Final_Site_JHU", 
                                              "Final_Site_MUSC", "Final_Site_NYU", "Final_Site_Rochester", 
                                              "Final_Site_UCSF", 'First_biop', 'Final_Chronicity', 
                                              "Final_ISN_[III][V]", "Final_ISN_[IV][V]", "Final_ISN_[V]"])


chronicityCNA_ISN_IVV = run_cna('Final_Activity', d, "SiteFBChron_ISN_IV_V", 
                                              ["Final_Site_Einstein", "Final_Site_JHU", 
                                              "Final_Site_MUSC", "Final_Site_NYU", "Final_Site_Rochester", 
                                              "Final_Site_UCSF", 'First_biop', 'Final_Chronicity', 
                                              "Final_ISN_[IV][V]"])

chronicityCNA_Pred = run_cna('Final_Activity', d, "SiteFBChron_Pred", 
                                              ["Final_Site_Einstein", "Final_Site_JHU", 
                                              "Final_Site_MUSC", "Final_Site_NYU", "Final_Site_Rochester", 
                                              "Final_Site_UCSF", 'First_biop', 'Final_Chronicity', 
                                              "Pred_use"])

chronicityCNA_Resp = run_cna('Final_Activity', d, "SiteFBChron_Resp", 
                                              ["Final_Site_Einstein", "Final_Site_JHU", 
                                              "Final_Site_MUSC", "Final_Site_NYU", "Final_Site_Rochester", 
                                              "Final_Site_UCSF", 'First_biop', 'Final_Chronicity', 
                                              "Responder_Status"])