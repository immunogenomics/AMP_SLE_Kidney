#!/bin/bash



cellbender remove-background --input /data/srlab2/public/srcollab/AMP_Phase_2/SLE/cellranger-5.0.1/GRCh38/AMPSLEkid_cells_0137/outs/raw_feature_bc_matrix.h5 --output /data/srlab2/public/srcollab/AMP_Phase_2/SLE/cellranger-5.0.1/GRCh38/AMPSLEkid_cells_0137/outs/raw_feature_bc_matrix_output.h5 --expected-cells 5000 --total-droplets-included 20000 --fpr 0.01 --epochs 150

