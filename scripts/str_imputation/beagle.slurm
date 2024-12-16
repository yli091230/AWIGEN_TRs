#!/bin/bash

#SBATCH --partition=ind-shared
#SBATCH --account=ddp268
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --time=24:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=32G
#SBATCH --export=ALL

module load cpu/0.15.4  gcc/10.2.0
module load anaconda3

set -e
chrom=$1
#log_dir="/expanse/protected/gymreklab-dbgap/mount/H3Africa/DS_I_Africa_project/from_ilifu/scripts/str_imputation/beagle_imputation_logs/"
impute_script="/expanse/protected/gymreklab-dbgap/mount/H3Africa/DS_I_Africa_project/from_ilifu/scripts/str_imputation/beagle_imputation.sh"
# create the folder if it is not exist
#mkdir -p ${log_dir}
time bash ${impute_script} ${chrom}
echo "All DONE!"
