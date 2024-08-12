#!/usr/bin/env python3

import argparse
import scrublet as scr
import scipy.io
import matplotlib.pyplot as plt
import numpy as np
import os
import pyreadr
import sys

parser = argparse.ArgumentParser()
parser.add_argument("run", type = str)
parser.add_argument("i", type = int)
args = parser.parse_args()
print(args.run)
print(args.i)

projected = np.load("/data/srlab2/jmears/jupyter/SLE_Phase2_Kidney_Analysis/cell-seq-output/python_data/cells_loaded.py")


counts_matrix = np.matrix.transpose(np.load("/data/srlab2/jmears/jupyter/SLE_Phase2_Kidney_Analysis/cell-seq-output/python_data/" + args.run + "_data_filtered_500nGene_1000nUMI_3pctnontargetMT.py"))
scrub = scr.Scrublet(counts_matrix, expected_doublet_rate=projected[args.i-1])


old_stdout = sys.stdout
sys.stdout = open('/data/srlab2/jmears/jupyter/SLE_Phase2_Kidney_Analysis/cell-seq-output/python_data/' + args.run + '_output.txt', 'a')
doublet_scores, predicted_doublets = scrub.scrub_doublets(min_counts=2,
                                                  min_cells=3,
                                                  min_gene_variability_pctl=90,
                                                  n_prin_comps=20)
sys.stdout.close()
sys.stdout = old_stdout


np.save("/data/srlab2/jmears/jupyter/SLE_Phase2_Kidney_Analysis/cell-seq-output/python_data/doublet_score_" + args.run, doublet_scores)

