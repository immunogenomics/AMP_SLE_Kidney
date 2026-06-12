#!/bin/bash
#SBATCH -J regressed_slingshot
#SBATCH -D /data/srlab/ssg34/SLE_kidney_v2/sugiarto_scripts/pseudotime_monoMacs/destiny/T_cells_v2/splitCD4and8
#SBATCH -o /data/srlab/ssg34/SLE_kidney_v2/sugiarto_scripts/pseudotime_monoMacs/destiny/T_cells_v2/splitCD4and8/destiny.out.txt
#SBATCH -e /data/srlab/ssg34/SLE_kidney_v2/sugiarto_scripts/pseudotime_monoMacs/destiny/T_cells_v2/splitCD4and8/destiny.err.txt
#SBATCH -p long,bigmem
#SBATCH --time=12:00:00
#SBATCH --mem=64G
#SBATCH -c 20
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=nsugiarto@bwh.harvard.edu

#------------------------ End of Header ------------------------#
source activate monocole3
Rscript runSlingshot.R