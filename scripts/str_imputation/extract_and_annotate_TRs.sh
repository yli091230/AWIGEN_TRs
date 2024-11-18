#!/bin/bash

: << 'comment'
need bcftools, tabix, TRTools
comment

chrom=$1

set -e

imputed_folder="/expanse/protected/gymreklab-dbgap/mount/H3Africa/DS_I_Africa_project/from_ilifu/results/imputed_TRs"
ref_file="/expanse/protected/gymreklab-dbgap/mount/yal084/Update_ensembleTR/reference_fixed_v4/ensembletr_refpanel_v4_chr${chrom}.vcf.gz"

echo "Extracting EnsembleTR IDs ... "
bcftools query -f '%ID\n' ${ref_file} | grep "EnsTR" > ${imputed_folder}/chr${chrom}_ensmble_ID.txt

echo "Using ID list to extract TRs ..."

tabix -p vcf ${imputed_folder}/chr${chrom}_imputed.vcf.gz 
bcftools view -i "ID=@${imputed_folder}/chr${chrom}_ensmble_ID.txt" ${imputed_folder}/chr${chrom}_imputed.vcf.gz -Oz -o ${imputed_folder}/chr${chrom}_imputed_TRs.vcf.gz

tabix -p vcf ${imputed_folder}/chr${chrom}_imputed_TRs.vcf.gz

echo "Add TR annotations ..."
annotaTR --vcf ${imputed_folder}/chr${chrom}_imputed_TRs.vcf.gz \
  --ref-panel ${ref_file} \
  --update-ref-alt \
  --vcftype hipstr \
  --dosages beagleap_norm \
  --out ${imputed_folder}/chr${chrom}_imputed_and_annotated_TRs \
  --outtype pgen vcf \
  --vcf-outtype z

echo "Successed!"
