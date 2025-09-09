#!/bin/bash/

folder=$(pwd)

mkdir -p log
declare -a types=("b_cell" "t_nk" "mono_dc")

for name in ${types[@]}; do
  mkdir -p cna_results/${name}
  echo $name
  cat << EOF | sbatch 
#!/bin/bash
#SBATCH -J ${name}_CNA
#SBATCH -D ${folder}
#SBATCH -o ${folder}/log/${name}_blood_CNA.out.txt
#SBATCH -e ${folder}/log/${name}_blood_CNA.err.txt
#SBATCH -p short
#SBATCH --time=3:00:00
#SBATCH --mem=8G
#SBATCH -c 8
#SBATCH --mail-type=end
#SBATCH --mail-user=$EMAIL

#------------------------ End of Header ------------------------#
source activate CTAPs

# python bloodCNA.py $name case_control
python bloodCNA.py $name chronicity
python bloodCNA.py $name activity

EOF
done