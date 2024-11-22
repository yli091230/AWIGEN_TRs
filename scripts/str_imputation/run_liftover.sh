#!/bin/bash

: << 'comment'
need to have pyliftover package installed,
need bcftools, tabix loaded
comment


LIFTOVER="/expanse/protected/gymreklab-dbgap/mount/H3Africa/DS_I_Africa_project/from_ilifu/scripts/str_imputation/liftover.py"
batch_size=1000
#input_file="/expanse/protected/gymreklab-dbgap/mount/H3Africa/DS_I_Africa_project/from_ilifu/results/converted_vcf/chr21_with_af.vcf.gz"
input_file="/expanse/protected/gymreklab-dbgap/mount/H3Africa/DS_I_Africa_project/from_ilifu/results/converted_vcf"
output_dir="/expanse/protected/gymreklab-dbgap/mount/H3Africa/DS_I_Africa_project/from_ilifu/results/converted_lifted_vcf"
include_SNPs="/expanse/protected/gymreklab-dbgap/mount/H3Africa/DS_I_Africa_project/from_ilifu/results/converted_vcf/SNPs_failed_QC.txt"

set -e
for chrom in $(seq 21 21); do
  tabix -p vcf ${input_file}/chr${chrom}_with_af.vcf.gz
  # remove SNPs where AF difference is too big between 1kg and target cohort
  bcftools view -R ${include_SNPs} ${input_file}/chr${chrom}_with_af.vcf.gz -Oz -o ${input_file}/chr${chrom}_with_af_filtered.vcf.gz
  mkdir -p ${output_dir}/chr${chrom}
  python3 ${LIFTOVER} ${input_file}/chr${chrom}_with_af_filtered.vcf.gz ${output_dir}/chr${chrom}
  bgzip ${output_dir}/chr${chrom}/chr${chrom}_with_af_filtered_hg38.vcf
  tabix -p vcf ${output_dir}/chr${chrom}/chr${chrom}_with_af_filtered_hg38.vcf.gz
  bgzip ${output_dir}/chr${chrom}/chr${chrom}_with_af_filtered_unlifted.vcf
  tabix -p vcf ${output_dir}/chr${chrom}/chr${chrom}_with_af_filtered_unlifted.vcf.gz
  # start split files
  # some of bcftools does not have the plugin function
  bcftools query -l ${output_dir}/chr${chrom}/chr${chrom}_with_af_filtered_hg38.vcf.gz | awk -v chrom=chr${chrom} -v group_size=${batch_size} 'BEGIN {FS=OFS="\t"} {print $1,"batch"int(NR / group_size)+1"_"chrom}' > ${output_dir}/chr${chrom}/chr${chrom}_sample_list.txt
  mkdir -p ${output_dir}/chr${chrom}/batch_files/
  awk -v outdir="${output_dir}/chr${chrom}/batch_files/" 'BEGIN {FS=OFS="\t"} {print $1 > outdir$2".txt"}' ${output_dir}/chr${chrom}/chr${chrom}_sample_list.txt
  mkdir -p ${output_dir}/chr${chrom}/split_by_samples
  for batch in ${output_dir}/chr${chrom}/batch_files/*.txt; do
    batch_name=$(basename ${batch} | cut -d "." -f 1)
    bcftools view -S ${batch} ${output_dir}/chr${chrom}/chr${chrom}_with_af_filtered_hg38.vcf.gz -Oz -o ${output_dir}/chr${chrom}/split_by_samples/${batch_name}.vcf.gz
  done
#  awk 'BEGIN {FS=OFS="\t"} {batch[$2] = batch[$2] ? batch[$2]","$1 : $1} END {for (b in batch) {print batch[b], b}}' ${output_dir}/chr${chrom}/chr${chrom}_sample_list.txt > ${output_dir}/chr${chrom}/chr${chrom}_sample_groups.txt
#	echo "get sample list"
#	bcftools plugin split ${output_dir}/chr${chrom}/chr${chrom}_with_af_filtered_hg38.vcf.gz -S ${output_dir}/chr${chrom}/chr${chrom}_sample_groups.txt -Oz -o ${output_dir}/chr${chrom}/split_by_samples/
  for f in ${output_dir}/chr${chrom}/split_by_samples/*.vcf.gz; do
    tabix -p vcf ${f}
  done
done


