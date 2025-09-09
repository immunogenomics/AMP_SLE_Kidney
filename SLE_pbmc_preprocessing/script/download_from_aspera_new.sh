module load aspera/default
export ASPERA_SCP_PASS=ImmportPass-40202

ascp -d -p -v -O 33001 -P 33001 --user=aspera.bapruzzese 'immport-upload15@aspera-immport.niaid.nih.gov:/AMP_RA_SLE.Phase2/RA\ PBMC\ RNA\ Seq/scRNAseq\ fastqs/10\ X\ Round\ 2\ RA-NU\December\ 2021' /data/srlab2/public/srcollab/AMP_Phase_2/SLE_pbmc/pilot_1121
