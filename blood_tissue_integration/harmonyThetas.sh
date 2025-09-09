#!/bin/bash/

folder=$(pwd)

mkdir -p log

for theta in {0..2}; do
  echo "Value: $theta"
  cat << EOF | sbatch 
#!/bin/bash
#SBATCH -J thetaLoop_$theta
#SBATCH -D ${folder}
#SBATCH -o ${folder}/log/thetaLoop_$theta.out.txt
#SBATCH -e ${folder}/log/thetaLoop_$theta.err.txt
#SBATCH -p short
#SBATCH --mem=16G
#SBATCH -c 8
#SBATCH --mail-type=end
#SBATCH --mail-user=$EMAIL

#------------------------ End of Header ------------------------#
# mkdir -p harmony_$theta

# source activate SLE

# Rscript startHarmony.R $theta

# source deactivate
source activate sc

# Rscript plot_UMAPs.R $theta
# Rscript split_UMAPs.R $theta

# Rscript nearestNeighbors.R $theta
Rscript bloodTissueScore.R $theta
EOF
done
