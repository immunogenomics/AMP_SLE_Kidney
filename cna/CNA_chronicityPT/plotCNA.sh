#!/bin/bash/

folder=$(pwd)

mkdir -p log
declare -a types=("t_nk" "myeloid" "GLOM" "DN" "LOH" "INTL")

for name in ${types[@]}; do
    if [ $name == "t_nk" ]; then
      nclus=5
      xPos=-3
      yPos=9
    elif [ $name == "myeloid" ]; then
      nclus=5
      xPos=-6.5
      yPos=-14
    elif [ $name == "GLOM" ]; then
      nclus=2
      xPos=-1
      yPos=-7.5
    elif [ $name == "PT" ]; then
      nclus=2
      xPos=-4.5
      yPos=-3.5
    elif [ $name == "DN" ]; then
      nclus=2
      xPos=1
      yPos=10
    elif [ $name == "LOH" ]; then
      nclus=3
      xPos=4.2
      yPos=-5.5
    elif [ $name == "INTL" ]; then
      nclus=3
      xPos=-2.5
      yPos=12
    fi
  mkdir -p cna_plots/${name}
  echo $name
  cat << EOF | sbatch 
#!/bin/bash
#SBATCH -J ${name}_CNA
#SBATCH -D ${folder}
#SBATCH -o ${folder}/log/${name}_plotCNA.out.txt
#SBATCH -e ${folder}/log/${name}_plotCNA.err.txt
#SBATCH -p short
#SBATCH --mem=8G
#SBATCH -c 8
#SBATCH --mail-type=end
#SBATCH --mail-user=$EMAIL

#------------------------ End of Header ------------------------#
source activate sc

Rscript plotCNA.R $name $nclus $xPos $yPos
EOF
done