#!/bin/bash

function do_mkfastq()
{
    local library=$1
    
cat << EOF | bsub 

#!/bin/bash
#BSUB -J ${library}
#BSUB -o /data/srlab2/public/srcollab/AMP_Phase_2/SLE_pbmc/pilot_1121/${library}/log/fastq_$library.out
#BSUB -e /data/srlab2/public/srcollab/AMP_Phase_2/SLE_pbmc/pilot_1121/${library}/log/fastq_$library.err
#BSUB -q big-multi
#BSUB -M 64000
#BSUB -n 16
#BSUB -N
#BSUB -R 'select[hname!=cn001]'
#BSUB -R 'select[hname!=cn002]'
#BSUB -R 'select[hname!=cn003]'
#BSUB -R 'select[hname!=cn004]'
#BSUB -R 'select[hname!=cn005]'

module load bcl2fastq2/patchelf-2.20.0
module load casava/1.8.3
cd /data/srlab2/public/srcollab/AMP_Phase_2/SLE_pbmc/pilot_1121/$library

TENX_IGNORE_DEPRECATED_OS=1 /data/srlab/cellranger/cellranger-6.1.1/cellranger mkfastq --id=FASTQS --run=/data/srlab2/public/srcollab/AMP_Phase_2/SLE_pbmc/pilot_1121/$library --csv=/data/srlab2/public/srcollab/AMP_Phase_2/SLE_pbmc/pilot_1121/$library/sample_sheet.csv --jobmode=local --localcores=16 --localmem=64

EOF
}

cat lsf_params_mkfast | while read library
do
#echo $library
do_mkfastq $library
done

