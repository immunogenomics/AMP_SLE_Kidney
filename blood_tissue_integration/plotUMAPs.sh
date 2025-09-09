#!/bin/bash/

folder=$(pwd)

mkdir -p log

declare -a types=("B" "Myeloid" "TNK")

for name in ${types[@]}; do
  mkdir -p $name
  echo $name
  theta=0
  # for theta in {0..2}; do
    mkdir -p $name/harmony_$theta
    echo "Value: $theta"
    # cat << EOF | sbatch 
#!/bin/bash
#SBATCH -J ${name}_${theta}
#SBATCH -D ${folder}
#SBATCH -o ${folder}/log/${name}_${theta}_UMAPs.out.txt
#SBATCH -e ${folder}/log/${name}_${theta}_UMAPs.err.txt
#SBATCH -p bigmem
#SBATCH --mem=32G
#SBATCH -c 8
#SBATCH --mail-type=end
#SBATCH --mail-user=$EMAIL

#------------------------ End of Header ------------------------#
# source activate SLE

# Rscript startHarmony.R $name $theta

# conda deactivate
# source activate sc

Rscript plot_UMAPs.R $name $theta
# Rscript split_UMAPs.R $name $theta

# Rscript nearestNeighbors.R $name $theta
Rscript bloodTissueScore.R $name $theta
# EOF
# done
done
