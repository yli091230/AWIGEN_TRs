#!/bin/bash

: <<'comment'
the --window determines the use of memory, default is 40. It requires at least 1.1x of --overlap parameters.
--overlap default value is 2.0
comment

chrom=$1
wid_size=40
#batch_size=1000
set -e
beagle="/expanse/protected/gymreklab-dbgap/mount/H3Africa/DS_I_Africa_project/from_ilifu/scripts/str_imputation/beagle.27May24.118.jar"

input_dir="/expanse/protected/gymreklab-dbgap/mount/H3Africa/DS_I_Africa_project/from_ilifu/results/converted_lifted_vcf/chr${chrom}/split_by_samples"
ensembleTR_ref="/expanse/protected/gymreklab-dbgap/mount/yal084/Update_ensembleTR/reference_fixed_v4"
genetic_map="/expanse/protected/gymreklab-dbgap/mount/H3Africa/DS_I_Africa_project/from_ilifu/additional_files/genetic_maps"
out_dir="/expanse/protected/gymreklab-dbgap/mount/H3Africa/DS_I_Africa_project/from_ilifu/results/imputed_TRs/chr${chrom}"
mkdir -p ${out_dir}/split_by_samples

# chr21_with_af_filtered_unlifted.vcf
# split samples by the batch_size

#bcftools query -l ${input_dir}/chr${chrom}_with_af_filtered_hg38.vcf | awk -v chrom=chr${chrom} -v group_size=${batch_size} 'BEGIN {FS=OFS="\t"} {print $1,"-","batch"int(NR / group_size)+1"_"chrom}' > ${out_dir}/chr${chrom}_sample_groups.txt
#bcftools plugin split ${input_dir}/chr${chrom}_with_af_filtered_hg38.vcf -G  
for f in ${input_dir}/*.vcf.gz; do
  file_name=$(basename ${f} | cut -d "." -f 1)
  # if current batch file finish imputed, skipping to next one
  if [[ -f ${out_dir}/split_by_samples/${file_name}_imputed_${wid_size}_window.vcf.gz.tbi ]]; then
    echo "${f} already imputed"
  else
    echo "Processing ${file_name}"
    # check if there any unfinished imputation file, if yes, delete the file before imputation
    if [[ -f ${out_dir}/split_by_samples/${file_name}_imputed_${wid_size}_window.vcf.gz ]]; then
      rm ${out_dir}/split_by_samples/${file_name}_imputed_${wid_size}_window.vcf.gz
    fi
  #  echo "Processing ${f}"
    java -Xmx24g -jar ${beagle} \
      gt=${f} \
      ref=${ensembleTR_ref}/ensembletr_refpanel_v4_chr${chrom}.bref3 \
      ap=true \
      map=${genetic_map}/plink.chr${chrom}.GRCh38_new.map \
      window=${wid_size} \
      out=${out_dir}/split_by_samples/${file_name}_imputed_${wid_size}_window
    tabix -p vcf ${out_dir}/split_by_samples/${file_name}_imputed_${wid_size}_window.vcf.gz
  fi
done
