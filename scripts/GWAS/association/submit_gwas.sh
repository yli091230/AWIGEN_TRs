

#variant_type=$1
read -p "Enter variant type (SNPs, TRs): " variant_type 
# check correct input
variant_list=("SNPs" "TRs")

if [[ ! "${variant_list[@]}" =~ "${variant_type}" ]];
then
  echo "Error: '$variant_type' is not in the list!" 
  echo "Please choose value from ${variant_list[@]}."
  exit 1
fi

gwas_script="/expanse/protected/gymreklab-dbgap/mount/H3Africa/DS_I_Africa_project/from_ilifu/scripts/GWAS/association/gwas.slurm"
log_folder="/expanse/protected/gymreklab-dbgap/mount/H3Africa/DS_I_Africa_project/from_ilifu/logs/gwas_association/hdl"

mkdir -p ${log_folder}
for chrom in $(seq 1 22); do
  sbatch --job-name=chr${chrom} --output=${log_folder}/chr${chrom}_gwas_$(date +%Y%m%d)_%j.out ${gwas_script} ${chrom} ${variant_type}
done
