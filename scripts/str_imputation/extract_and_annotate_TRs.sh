#!/bin/bash

: << 'comment'
need bcftools, tabix, TRTools
e.g. to extract only TR from chr22 run:
./extract_and_annotate_TRs.sh 22 true
comment

chrom=$1
str_only=$2

# Check if str_only is a valid value: "true" or "false"
if [[ "$str_only" != "true" && "$str_only" != "false" ]]; then
    echo "A valid value for the second argument needs to be: true or false"
fi

window_size=40
set -e

imputed_folder="/expanse/protected/gymreklab-dbgap/mount/H3Africa/DS_I_Africa_project/from_ilifu/results/imputed_TRs/chr${chrom}"
ref_file="/expanse/protected/gymreklab-dbgap/mount/yal084/Update_ensembleTR/reference_fixed_v4/ensembletr_refpanel_v4_chr${chrom}.vcf.gz"

mkdir -p ${imputed_folder}
if [[ "$str_only" == "true" ]]; then

  echo "Extracting EnsembleTR IDs ... "
  bcftools query -f '%ID\n' ${ref_file} | grep "EnsTR" > ${imputed_folder}/chr${chrom}_ensmble_ID.txt

  echo "Using ID list to extract TRs ..."

  #tabix -p vcf ${imputed_folder}/chr${chrom}_imputed.vcf.gz
  for vcf in ${imputed_folder}/split_by_samples/batch*_chr${chrom}_imputed_${window_size}_window.vcf.gz; do 
    vcf_name=$(basename ${vcf} | cut -d "." -f 1)
    bcftools view -i "ID=@${imputed_folder}/chr${chrom}_ensmble_ID.txt" ${vcf} -Oz -o ${imputed_folder}/split_by_samples/${vcf_name}_TRs.vcf.gz
  #bcftools view -i "ID=@${imputed_folder}/chr${chrom}_ensmble_ID.txt" ${imputed_folder}/chr${chrom}_imputed.vcf.gz -Oz -o ${imputed_folder}/chr${chrom}_imputed_TRs.vcf.gz
    tabix -p vcf ${imputed_folder}/split_by_samples/${vcf_name}_TRs.vcf.gz
  done
  echo "finish extracting TRs for chr${chrom}, start merge TRs ..."
  echo ${imputed_folder}/split_by_samples/*_TRs.vcf.gz

  # merge all samples together
  bcftools merge ${imputed_folder}/split_by_samples/*_TRs.vcf.gz -Oz -o ${imputed_folder}/chr${chrom}_imputed_TRs.vcf.gz
  bcftools index -t -f ${imputed_folder}/chr${chrom}_imputed_TRs.vcf.gz
  merged_vcf="${imputed_folder}/chr${chrom}_imputed_TRs.vcf.gz"
  out_prefix="${imputed_folder}/chr${chrom}_imputed_and_annotated_TRs"
else
  echo "Use all variants, starting combining..."
  ls ${imputed_folder}/split_by_samples/batch*_chr${chrom}_imputed_${window_size}_window.vcf.gz
  bcftools merge ${imputed_folder}/split_by_samples/batch*_chr${chrom}_imputed_${window_size}_window.vcf.gz -Oz -o ${imputed_folder}/chr${chrom}_imputed_merged_all_variants.vcf.gz
  bcftools index -t -f ${imputed_folder}/chr${chrom}_imputed_merged_all_variants.vcf.gz
  merged_vcf="${imputed_folder}/chr${chrom}_imputed_merged_all_variants.vcf.gz"
  out_prefix="${imputed_folder}/chr${chrom}_imputed_and_annotated_all_variants"
fi

# may need to check the order of samples to make sure it consistent across different files
echo "Add TR annotations ..."
annotaTR --vcf ${merged_vcf} \
  --ref-panel ${ref_file} \
  --update-ref-alt \
  --vcftype hipstr \
  --dosages beagleap_norm \
  --out ${out_prefix} \
  --outtype pgen vcf \
  --vcf-outtype z

echo "Successed!"
