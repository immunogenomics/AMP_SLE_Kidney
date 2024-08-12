#!/bin/bash

function do_count()
{
    local library=$1
        local sample=$2
    local transcriptome=$3
	local features=$4
    local libraries=$5
        local version=$6
    local genome=$7

cat << EOF | bsub 

#!/bin/bash
#BSUB -J ${sample}
#BSUB -o /data/srlab2/public/srcollab/AMP_Phase_2/SLE_pbmc/pilot_1121/${library}/log/${sample}.out
#BSUB -e /data/srlab2/public/srcollab/AMP_Phase_2/SLE_pbmc/pilot_1121/${library}/log/${sample}.err
#BSUB -q big-multi
#BSUB -M 32000
#BSUB -n 8
#BSUB -N
#BSUB -R 'select[hname!=cn001]'
#BSUB -R 'select[hname!=cn002]'
#BSUB -R 'select[hname!=cn003]'
#BSUB -R 'select[hname!=cn004]'
#BSUB -R 'select[hname!=cn005]'

cd /data/srlab2/public/srcollab/AMP_Phase_2/SLE_pbmc/pilot_1121/$library/$version/$genome


TENX_IGNORE_DEPRECATED_OS=1 /data/srlab/cellranger/cellranger-6.1.1/cellranger count --id=$sample --libraries=$libraries --feature-ref=$features --transcriptome=$transcriptome --jobmode=local --localcores=16 --localmem=64

EOF
}

cat lsf_params_count_features_1121_NW | while read library sample transcriptome features libraries version genome
do
do_count $library $sample $transcriptome $features $libraries $version $genome
done
