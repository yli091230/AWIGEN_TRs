

beagle_script="/expanse/protected/gymreklab-dbgap/mount/H3Africa/DS_I_Africa_project/from_ilifu/scripts/str_imputation/beagle.slurm"
log_folder="/expanse/protected/gymreklab-dbgap/mount/H3Africa/DS_I_Africa_project/from_ilifu/logs/beagle_imputation"

for chrom in $(seq 1 22); do
  mkdir -p ${log_folder}/chr${chrom}
  sbatch --job-name=chr${chrom} --output=${log_folder}/chr${chrom}/chr${chrom}_imputation_$(date +%Y%m%d)_%j.out ${beagle_script} ${chrom}

done
