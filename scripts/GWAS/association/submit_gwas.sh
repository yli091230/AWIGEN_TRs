

gwas_script="/expanse/protected/gymreklab-dbgap/mount/H3Africa/DS_I_Africa_project/from_ilifu/scripts/GWAS/association/gwas.slurm"
log_folder="/expanse/protected/gymreklab-dbgap/mount/H3Africa/DS_I_Africa_project/from_ilifu/logs/gwas_association/ldl_c"

mkdir -p ${log_folder}
for chrom in $(seq 19 19); do
  sbatch --job-name=chr${chrom} --output=${log_folder}/chr${chrom}_gwas_$(date +%Y%m%d)_%j.out ${gwas_script} ${chrom}
done
