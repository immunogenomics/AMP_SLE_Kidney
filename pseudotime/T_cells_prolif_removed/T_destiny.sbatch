#!/bin/bash
#SBATCH -J noProlif_pseudotime
#SBATCH -D /data/srlab/ssg34/SLE_kidney_v2/sugiarto_scripts/pseudotime_monoMacs/destiny/T_cells_noProlif
#SBATCH -o /data/srlab/ssg34/SLE_kidney_v2/sugiarto_scripts/pseudotime_monoMacs/destiny/T_cells_noProlif/logs/Tcell.out.txt
#SBATCH -e /data/srlab/ssg34/SLE_kidney_v2/sugiarto_scripts/pseudotime_monoMacs/destiny/T_cells_noProlif/logs/Tcell.err.txt
#SBATCH -p long,bigmem
#SBATCH --time=12:00:00
#SBATCH --mem=64G
#SBATCH -c 20
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=nsugiarto@bwh.harvard.edu

#------------------------ End of Header ------------------------#
source activate destiny

Rscript T_destiny.R

sbatch T_slingshot.sbatch