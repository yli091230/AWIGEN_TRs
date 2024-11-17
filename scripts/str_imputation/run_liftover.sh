#!/bin/bash

LIFTOVER="/expanse/protected/gymreklab-dbgap/mount/H3Africa/DS_I_Africa_project/from_ilifu/scripts/str_imputation/liftover.py"
#input_file="/expanse/protected/gymreklab-dbgap/mount/H3Africa/DS_I_Africa_project/from_ilifu/results/converted_vcf/chr21_with_af.vcf.gz"
input_file="/expanse/protected/gymreklab-dbgap/mount/H3Africa/DS_I_Africa_project/from_ilifu/results/converted_vcf"
output_dir="/expanse/protected/gymreklab-dbgap/mount/H3Africa/DS_I_Africa_project/from_ilifu/results/converted_lifted_vcf"

for chrom in $(seq 1 22); do
  python3 ${LIFTOVER} ${input_file}/chr${chrom}_with_af.vcf.gz ${output_dir}
done


