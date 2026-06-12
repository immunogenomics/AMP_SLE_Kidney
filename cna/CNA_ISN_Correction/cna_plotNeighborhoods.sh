#!/bin/bash/

folder=$(pwd)

mkdir -p log
declare -a types=("b_plasma" "t_nk" "myeloid" "GLOM" "PT" "DN" "LOH" "INTL")

for name in ${types[@]}; do
  if [ $name == "t_nk" ]; then
    meta='/data/srlab/ssg34/SLE_kidney_v2/data/qcd/T_NK_clusterQCd_meta_harmonizedPCUMAPclusters_annotations02062024.rds'
  elif [ $name == "b_plasma" ]; then
    meta='/data/srlab/ssg34/SLE_kidney_v2/data/qcd/BP_clusterQCd_cellstate_meta_annotations_09232022.rds'
  elif [ $name == "myeloid" ]; then
    meta='/data/srlab/ssg34/SLE_kidney_v2/data/qcd/Myeloid_clusterQCd_meta_harmonizedPCUMAPCellStateClusters_10042022.rds'
  elif [ $name == "GLOM" ]; then
    meta='/data/srlab/ssg34/SLE_kidney_v2/data/tissue/glom_meta_qcd_harmony_umap_clusternames_11302023.rds'
  elif [ $name == "PT" ]; then
    meta='/data/srlab/ssg34/SLE_kidney_v2/data/tissue/pt_meta_qcd_harmony_umap_clusternames_11302023.rds'
  elif [ $name == "DN" ]; then
    meta='/data/srlab/ssg34/SLE_kidney_v2/data/tissue/dn_meta_qcd_harmony_umap_clusternames_11302023.rds'
  elif [ $name == "LOH" ]; then
    meta='/data/srlab/ssg34/SLE_kidney_v2/data/tissue/loh_meta_qcd_harmony_umap_clusternames_11302023.rds'
  elif [ $name == "INTL" ]; then
    meta='/data/srlab/ssg34/SLE_kidney_v2/data/tissue/intl_meta_qcd_harmony_umap_clusternames_11302023.rds'
  fi
  echo $meta
  for file in cna_results/${name}_tissue/*_ncorr.csv; do
    mkdir -p cna_plots/${name}_tissue
    echo $file
    cat << EOF | sbatch 
#!/bin/bash
#SBATCH -J ${name}_CNA
#SBATCH -D ${folder}
#SBATCH -o ${folder}/log/${name}_tissue_$(basename $file)_plotEm.out.txt
#SBATCH -e ${folder}/log/${name}_tissue_$(basename $file)_plotEm.err.txt
#SBATCH -p short
#SBATCH --mem=8G
#SBATCH -c 8
#SBATCH --mail-type=end
#SBATCH --mail-user=$EMAIL

#------------------------ End of Header ------------------------#
source activate sc

Rscript plotNeighborhoods.R $name tissue $file $meta
EOF
  done
done

declare -a types=("b_cell" "t_nk" "mono_dc")

for name in ${types[@]}; do
  if [ $name == "t_nk" ]; then
    meta='/data/srlab/ssg34/SLE_kidney_v2/data/pbmc/SLE_t_nk_pbmc_metadata_filtered_V0_cases_controls_07182024.rds'
  elif [ $name == "b_cell" ]; then
    meta='/data/srlab/ssg34/SLE_kidney_v2/data/pbmc/SLE_b_cell_pbmc_metadata_filtered_V0_cases_controls_07182024.rds'
  elif [ $name == "mono_dc" ]; then
    meta='/data/srlab/ssg34/SLE_kidney_v2/data/pbmc/SLE_mono_dc_pbmc_metadata_filtered_V0_cases_controls_07182024.rds'
  fi
  echo $name
  for file in cna_results/${name}_pbmc/*_ncorr.csv; do
    mkdir -p cna_plots/${name}_pbmc
    echo $file
    cat << EOF | sbatch 
#!/bin/bash
#SBATCH -J ${name}_CNA
#SBATCH -D ${folder}
#SBATCH -o ${folder}/log/${name}_pbmc_$(basename $file)_plotEm.out.txt
#SBATCH -e ${folder}/log/${name}_pbmc_$(basename $file)_plotEm.err.txt
#SBATCH -p short
#SBATCH --mem=8G
#SBATCH -c 8
#SBATCH --mail-type=end
#SBATCH --mail-user=$EMAIL

#------------------------ End of Header ------------------------#
source activate sc

Rscript plotNeighborhoods.R $name pbmc $file $meta
EOF
  done
done