#!/bin/bash

: <<'comment'
the --window determines the use of memory, default is 40. It requires at least 1.1x of --overlap parameters.
--overlap default value is 2.0
comment

chrom=$1
set -e
beagle="/expanse/protected/gymreklab-dbgap/mount/H3Africa/DS_I_Africa_project/from_ilifu/scripts/str_imputation/beagle.27May24.118.jar"

input_dir="/expanse/protected/gymreklab-dbgap/mount/H3Africa/DS_I_Africa_project/from_ilifu/results/converted_lifted_vcf"
ensembleTR_ref="/expanse/protected/gymreklab-dbgap/mount/yal084/Update_ensembleTR/reference_fixed_v4"
genetic_map="/expanse/protected/gymreklab-dbgap/mount/H3Africa/DS_I_Africa_project/from_ilifu/additional_files/genetic_maps"
out_dir="/expanse/protected/gymreklab-dbgap/mount/H3Africa/DS_I_Africa_project/from_ilifu/results/imputed_TRs"

java -Xmx24g -jar ${beagle} \
  gt=${input_dir}/chr${chrom}_with_af_hg38.vcf \
  ref=${ensembleTR_ref}/ensembletr_refpanel_v4_chr${chrom}.bref3 \
  ap=true \
  map=${genetic_map}/plink.chr${chrom}.GRCh38_new.map \
  window=5 \
  out=${out_dir}/chr${chrom}_imputed 
#tabix -p vcf

