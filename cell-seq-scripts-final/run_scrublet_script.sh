#!/bin/bash

function do_job()
{
    local run=$1
       local i=$2

cat << EOF | bsub 

#!/bin/bash
#BSUB -J save_raw_${run}
#BSUB -o /data/srlab2/jmears/jupyter/SLE_Phase2_Kidney_Analysis/cell-seq-scripts-final/log/scrublet_${run}.out
#BSUB -e /data/srlab2/jmears/jupyter/SLE_Phase2_Kidney_Analysis/cell-seq-scripts-final/log/scrublet_${run}.err
#BSUB -q normal
#BSUB -R 'select[hname!=cn001]'
#BSUB -R 'select[hname!=cn002]'
#BSUB -R 'select[hname!=cn003]'
#BSUB -R 'select[hname!=cn004]'
#BSUB -R 'select[hname!=cn005]'

./scrublet_script.py $run $i


EOF
}

cat list_of_cell_samples.txt | while read run i
do
do_job $run $i
done
