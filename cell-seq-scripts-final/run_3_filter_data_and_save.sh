#!/bin/bash

function do_job()
{
    local run=$1

cat << EOF | bsub 

#!/bin/bash
#BSUB -J save_raw_${run}
#BSUB -o /data/srlab2/jmears/jupyter/SLE_Phase2_Kidney_Analysis/cell-seq-scripts-final/log/filter_${run}.out
#BSUB -e /data/srlab2/jmears/jupyter/SLE_Phase2_Kidney_Analysis/cell-seq-scripts-final/log/filter_${run}.err
#BSUB -q bigmem
#BSUB -R 'select[hname!=cn001]'
#BSUB -R 'select[hname!=cn002]'
#BSUB -R 'select[hname!=cn003]'
#BSUB -R 'select[hname!=cn004]'
#BSUB -R 'select[hname!=cn005]'

./3_filter_data_and_save.R ${run}


EOF
}

cat list_of_cell_samples.txt | while read run num
do
do_job $run
done
