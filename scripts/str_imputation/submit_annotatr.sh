

beagle_script="/expanse/protected/gymreklab-dbgap/mount/H3Africa/DS_I_Africa_project/from_ilifu/scripts/str_imputation/annotatr.slurm"
log_folder="/expanse/protected/gymreklab-dbgap/mount/H3Africa/DS_I_Africa_project/from_ilifu/logs/merge_and_annotate_TRs/"
str_only="true"

for chrom in $(seq 21 21); do
  mkdir -p ${log_folder}/chr${chrom}
  sbatch --job-name=chr${chrom} --output=${log_folder}/chr${chrom}/chr${chrom}_merge_and_annotation_$(date +%Y%m%d)_%j.out ${beagle_script} ${chrom} ${str_only}

done
