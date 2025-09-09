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
type=sys.argv[2]

harmony = pd.read_csv(('/data/srlab/ssg34/SLE_kidney_v2/sugiarto_scripts/CNA_bloodOverlap/cna_inputs/' + name + '/' + type + '/harmony.csv'))
meta = pd.read_csv(('/data/srlab/ssg34/SLE_kidney_v2/sugiarto_scripts/CNA_bloodOverlap/cna_inputs/' + name + '/' + type + '/meta.csv'))
umap = pd.read_csv(('/data/srlab/ssg34/SLE_kidney_v2/sugiarto_scripts/CNA_bloodOverlap/cna_inputs/' + name + '/' + type + '/umap.csv'))

if type == "chronicity":
    var = "Final_Chronicity"
    covars = ["First_biop", 
                 "Final_Site_Einstein", "Final_Site_JHU", 
                 "Final_Site_MUSC", "Final_Site_NYU", "Final_Site_Rochester", 
                 "Final_Site_UCSF", "Sex", "Age"]
elif type == "activity":
    var = "Final_Activity"
    covars = ["Final_Site_Einstein", "Final_Site_JHU", 
                 "Final_Site_MUSC", "Final_Site_NYU", "Final_Site_Rochester", 
                 "Final_Site_UCSF", 'First_biop', 'Final_Chronicity', "Sex", "Age"]
elif type == "case_control":
    var = "Type"
    covars = ['Sex_numeric', 'Age']


d = mad.MultiAnnData(X=harmony, obs=meta, sampleid="Unified_Visit")

d.obs_to_sample(([var] + covars))

umap.index = d.obs.index
d.obsm['X_umap'] = umap

cna.pp.knn(d)

with open((name + "_" + type + "_pValues.csv"), 'w') as file:
    # Append the line of text
    file.write("name,pValue\n")

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
    with open((name + "_" + type + "_pValues.csv"), 'a') as file:
        # Append the line of text
        file.write(fileName + "," + str(cna_res.p) + '\n')

    np.savetxt("cna_results/" + name + "/" + variable + "_" + fileName + "_ncorr.csv", 
                cna_res.ncorrs, delimiter=",")
    np.savetxt("cna_results/" + name + "/" + variable + "_" + fileName + "_fdrs.csv", 
                cna_res.fdrs, delimiter=",")

run_cna(var, d, "none")
run_cna(var, d, "tissueCorrection", covars[:-2])
run_cna(var, d, "sexAge", ["Sex", "Age"])