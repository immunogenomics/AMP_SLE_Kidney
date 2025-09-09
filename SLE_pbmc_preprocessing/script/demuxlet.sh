#!/bin/bash

#BSUB -o /data/srlab2/public/srcollab/AMP_Phase_2/SLE_pbmc/genotyping/BROAD/pilot/%J.out
#BSUB -e /data/srlab2/public/srcollab/AMP_Phase_2/SLE_pbmc/genotyping/BROAD/pilot/%J.out
#BSUB -J demux_broad
#BSUB -q big-multi
#BSUB -M 64000
#BSUB -n 16

./run_demuxlet.sh
