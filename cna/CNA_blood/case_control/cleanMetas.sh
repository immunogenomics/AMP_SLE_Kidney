#!/bin/bash/

folder=$(pwd)

mkdir -p log
declare -a types=("b_cell" "t_nk" "mono_dc")

for name in ${types[@]}; do

  mkdir -p new_metas/$name

  Rscript cleanBloodMeta.R $name case_control
done