#!/bin/bash


PLINK2="/expanse/protected/gymreklab-dbgap/mount/H3Africa/DS_I_Africa_project/from_ilifu/scripts/plinkTovcf/plink_2.0/plink2"

all_snps="/expanse/protected/gymreklab-dbgap/mount/H3Africa/DS_I_Africa_project/from_ilifu/awigen_dataset/plink_file/awigen-qc"

filtered_snps="/expanse/protected/gymreklab-dbgap/mount/H3Africa/DS_I_Africa_project/from_ilifu/results/GWAS/PCs/filtering_SNPs/biallelic"

kinship_folder="/expanse/protected/gymreklab-dbgap/mount/H3Africa/DS_I_Africa_project/from_ilifu/results/GWAS/kinship/"

mkdir -p ${kinship_folder}

#echo "converting to binary"
#${PLINK2} --file ${input_file} \
#  --make-bed \
#  --out ${input_file}_binary

#echo "calculate kinship"

echo "   Using the filtered SNPs as input"
${PLINK2} --bfile ${filtered_snps} \
  --make-king-table \
  --out ${kinship_folder}/filtered_SNPs_kinship

# echo "   Using all SNPs as input"
# ${PLINK2} --bfile ${all_snps} \
#   --make-king-table \
#   --out ${kinship_folder}/all_SNPs_kinship

