#!/bin/bash

chrom=$1
variant_type="SNPs"
PLINK2="/expanse/protected/gymreklab-dbgap/mount/H3Africa/DS_I_Africa_project/from_ilifu/scripts/plinkTovcf/plink_2.0/plink2"

#pgen_file="/expanse/protected/gymreklab-dbgap/mount/H3Africa/DS_I_Africa_project/from_ilifu/results/imputed_TRs/chr${chrom}/chr${chrom}_imputed_and_annotated_TRs"
pgen_file="/expanse/protected/gymreklab-dbgap/mount/H3Africa/DS_I_Africa_project/from_ilifu/results/imputed_TRs/chr${chrom}/chr${chrom}_imputed_filtered_${variant_type}"

pheno_file="/expanse/protected/gymreklab-dbgap/mount/H3Africa/DS_I_Africa_project/from_ilifu/results/GWAS/regression_files/friedewald_ldl_c/friedewald_ldl_pheno.tsv"
cov_file="/expanse/protected/gymreklab-dbgap/mount/H3Africa/DS_I_Africa_project/from_ilifu/results/GWAS/regression_files/friedewald_ldl_c/ldlc_covs.tsv"
kin_table="/expanse/protected/gymreklab-dbgap/mount/H3Africa/DS_I_Africa_project/from_ilifu/results/GWAS/kinship/filtered_${variant_type}_kinship.kin0"

out_file="/expanse/protected/gymreklab-dbgap/mount/H3Africa/DS_I_Africa_project/from_ilifu/results/GWAS/associations/friedewald_ldl_c/"

mkdir -p ${out_file}

${PLINK2} --pfile ${pgen_file} \
  --pheno "iid-only" ${pheno_file} \
  --glm hide-covar \
  --covar "iid-only" ${cov_file} \
  --covar-variance-standardize \
  --mac 20 \
  --ci 0.95 \
  --king-cutoff-table ${kin_table} 0.177 \
  --out ${out_file}/merged_chr${chrom}_${variant_type}_gwas 
