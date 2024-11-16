#!/bin/bash

: <<'comment'
  This script convert the plink files to VCF:
  1. Use the hg19 reference to set REF/ALT alleles
  2. Only keep snps with {A,C,G,T,a,c,g,t,missing}
  3. calculate the allele frequency

  Required module:
  bcftools
comment

#chrom=21
PLINK="/expanse/protected/gymreklab-dbgap/mount/H3Africa/DS_I_Africa_project/from_ilifu/scripts/plinkTovcf/plink_1.9/plink"
PLINK2="/expanse/protected/gymreklab-dbgap/mount/H3Africa/DS_I_Africa_project/from_ilifu/scripts/plinkTovcf/plink_2.0/plink2"
HG19="/expanse/protected/gymreklab-dbgap/mount/H3Africa/DS_I_Africa_project/from_ilifu/additional_files/reference_genome/human_g1k_v37.fasta"
plink_file="/expanse/protected/gymreklab-dbgap/mount/H3Africa/DS_I_Africa_project/from_ilifu/plink_file/awigen-qc"
output_dir="/expanse/protected/gymreklab-dbgap/mount/H3Africa/DS_I_Africa_project/from_ilifu/results/converted_vcf"

set -e
for chrom in $(seq 1 22); do 
  echo "Processing chrom${chrom}"
  ${PLINK2} \
    --bfile ${plink_file} \
    --chr ${chrom} \
    --recode vcf id-paste=iid bgz \
    --snps-only 'just-acgt' \
    --fa ${HG19} \
    --ref-from-fa \
    --out ${output_dir}/chr${chrom}
  ls ${output_dir}/chr${chrom}.vcf.gz
  bcftools +fill-tags ${output_dir}/chr${chrom}.vcf.gz -Oz -o ${output_dir}/chr${chrom}_with_af.vcf.gz -- -t AF
  ls ${output_dir}/chr${chrom}_with_af.vcf.gz
  rm ${output_dir}/chr${chrom}.vcf.gz
  bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t%AF\n' ${output_dir}/chr${chrom}_with_af.vcf.gz > ${output_dir}/chr${chrom}_allele_frequencies.txt
# bcftools +fill-tags input.vcf -- -t CALLRATE
done
