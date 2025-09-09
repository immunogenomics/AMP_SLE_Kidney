#!/bin/bash/

folder=$(pwd)

mkdir -p log
declare -a types=("B" "TNK" "Myeloid")


for name in ${types[@]}; do

  if [ $name == "TNK" ]; then
    tissue_norm='/data/srlab/ssg34/SLE_kidney_v2/data/qcd/T_NK_clusterQCd_norm_08112024.rds'
    tissue_meta='/data/srlab/ssg34/SLE_kidney_v2/data/qcd/T_NK_clusterQCd_meta_harmonizedPCUMAPclusters_annotations02062024.rds'
    pbmc_norm='/data/srlab/ssg34/SLE_pbmc_analysis/data/gene_expression/t_nk_qcd_seurat_sc_analysis_08252023.rds'
    pbmc_meta='/data/srlab/ssg34/SLE_kidney_v2/data/pbmc/SLE_t_nk_pbmc_metadata_filtered_V0_cases_controls_07182024.rds'
    ref='/data/srlab/ssg34/SLE_kidney_v2/data/symphony_materials/t_nk_kidney_symphony_reference_07162024.rds'
  elif [ $name == "B" ]; then
    tissue_norm='/data/srlab/ssg34/SLE_kidney_v2/data/qcd/BP_clusterQCd_norm_09232022.rds'
    tissue_meta='/data/srlab/ssg34/SLE_kidney_v2/data/qcd/BP_clusterQCd_cellstate_meta_annotations_09232022.rds'
    pbmc_norm='/data/srlab/ssg34/SLE_pbmc_analysis/data/gene_expression/b_cell_qcd_seurat_sc_analysis_08252023.rds'
    pbmc_meta='/data/srlab/ssg34/SLE_kidney_v2/data/pbmc/SLE_b_cell_pbmc_metadata_filtered_V0_cases_controls_07182024.rds'
    ref='/data/srlab/ssg34/SLE_kidney_v2/data/symphony_materials/bp_kidney_symphony_reference_07162024.rds'
  elif [ $name == "Myeloid" ]; then
    tissue_norm='/data/srlab/ssg34/SLE_kidney_v2/data/qcd/Myeloid_clusterQCd_norm_10042022.rds'
    tissue_meta='/data/srlab/ssg34/SLE_kidney_v2/data/qcd/Myeloid_clusterQCd_meta_harmonizedPCUMAPCellStateClusters_10042022.rds'
    pbmc_norm='/data/srlab/ssg34/SLE_pbmc_analysis/data/gene_expression/mono_dc_qcd_seurat_sc_analysis_08252023.rds'
    pbmc_meta='/data/srlab/ssg34/SLE_kidney_v2/data/pbmc/SLE_mono_dc_pbmc_metadata_filtered_V0_cases_controls_07182024.rds'
    ref='/data/srlab/ssg34/SLE_kidney_v2/data/symphony_materials/myeloid_kidney_symphony_reference_07162024.rds'
  fi
  mkdir -p PCs/
  echo $name
  # cat << EOF | sbatch 
#!/bin/bash
#SBATCH -J ${name}_PCA
#SBATCH -D ${folder}
#SBATCH -o ${folder}/log/${name}_PCA.out.txt
#SBATCH -e ${folder}/log/${name}_PCA.err.txt
#SBATCH -p long
#SBATCH --mem=64G
#SBATCH -c 8
#SBATCH --mail-type=end
#SBATCH --mail-user=$EMAIL

#------------------------ End of Header ------------------------#
# source activate sc

Rscript weightedPCA.R $tissue_norm $tissue_meta $pbmc_norm $pbmc_meta $ref $name
# EOF
done
