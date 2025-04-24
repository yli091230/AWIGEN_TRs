#!/bin/bash

#SBATCH --partition=ind-shared
#SBATCH --account=ddp268
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --time=8:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=2G
#SBATCH --export=ALL

module load cpu/0.15.4  gcc/10.2.0
module load anaconda3
module load tabix
module load bcftools
#eval "$(conda shell.bash hook)"
#conda activate trtools

set -e
chrom=$1
variant_type=$2
gwas_script="/expanse/protected/gymreklab-dbgap/mount/H3Africa/DS_I_Africa_project/from_ilifu/scripts/GWAS/association/run_associ.sh"
time bash ${gwas_script} ${chrom} ${variant_type}
echo "All DONE!"
