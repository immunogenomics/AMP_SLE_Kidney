#!/bin/bash

function do_mkfastq()
{
    local library=$1
    
cat << EOF | bsub 

#!/bin/bash
#BSUB -J ${library}
#BSUB -o /data/srlab2/public/srcollab/AMP_Phase_2/SLE_pbmc/ADT/high_yield_Nextseq/${library}/log/fastq_$library.out
#BSUB -e /data/srlab2/public/srcollab/AMP_Phase_2/SLE_pbmc/ADT/high_yield_Nextseq/${library}/log/fastq_$library.err
#BSUB -q big-multi
#BSUB -M 32000
#BSUB -n 8
#BSUB -N
#BSUB -R 'select[hname!=cn001]'
#BSUB -R 'select[hname!=cn002]'
#BSUB -R 'select[hname!=cn003]'
#BSUB -R 'select[hname!=cn004]'
#BSUB -R 'select[hname!=cn005]'

module load bcl2fastq2/2.19.1
module load casava/1.8.3
cd /data/srlab2/public/srcollab/AMP_Phase_2/SLE_pbmc/ADT/high_yield_Nextseq/$library
 

/data/srlab/cellranger/cellranger-6.0.1/cellranger mkfastq --id=FASTQS_NEW --run=/data/srlab2/public/srcollab/AMP_Phase_2/SLE_pbmc/ADT/high_yield_Nextseq/$library --csv=/data/srlab2/public/srcollab/AMP_Phase_2/SLE_pbmc/ADT/high_yield_Nextseq/$library/thesamplesheet.csv --lanes=1,2,3,4 --jobmode=local --localcores=16 --localmem=64

EOF
}

cat lsf_params_mkfast_cite | while read library
do
#echo $library
do_mkfastq $library
done

