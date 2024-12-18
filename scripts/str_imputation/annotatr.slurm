#!/bin/bash

#SBATCH --partition=ind-shared
#SBATCH --account=ddp268
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --time=12:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=2G
#SBATCH --export=ALL

module load cpu/0.15.4  gcc/10.2.0
module load anaconda3
module load tabix
module load bcftools

conda activate trtools

set -e
chrom=$1
str_only="true"
anno_script="/expanse/protected/gymreklab-dbgap/mount/H3Africa/DS_I_Africa_project/from_ilifu/scripts/str_imputation/extract_and_annotate_TRs.sh"
time bash ${anno_script} ${chrom} ${str_only}
echo "All DONE!"
