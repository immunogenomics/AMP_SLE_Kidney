#/bin/bash 
for i in channel_1 channel_2 channel_3

do

conda activate demuxlet
mkdir /data/srlab2/public/srcollab/AMP_Phase_2/SLE_pbmc/pilot_1121/BroadRunin2/genotyping/$i

cmd="demuxlet --sam /data/srlab2/public/srcollab/AMP_Phase_2/SLE_pbmc/pilot_1121/BroadRunin2/cellranger-6.1.1/GRCh38/$i/outs/possorted_genome_bam.bam --tag-group CB --tag-UMI UB --field GT --vcf /data/srlab2/public/srcollab/AMP_Phase_2/SLE_pbmc/vcf/BROAD_Pilot/filtered_ordered_LUPUS_PBMC.vcf.gz --out /data/srlab2/public/srcollab/AMP_Phase_2/SLE_pbmc/pilot_1121/BroadRunin2/genotyping/$i"

bsubhosts -M 30000 -n 10 $cmd
done
