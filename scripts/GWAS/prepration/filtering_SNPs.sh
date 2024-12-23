
PLINK2="/expanse/protected/gymreklab-dbgap/mount/H3Africa/DS_I_Africa_project/from_ilifu/scripts/plinkTovcf/plink_2.0/plink2"
MINMAF=0.05
out_prefix=/expanse/protected/gymreklab-dbgap/mount/H3Africa/DS_I_Africa_project/from_ilifu/results/GWAS/PCs

bfile="/expanse/protected/gymreklab-dbgap/mount/H3Africa/DS_I_Africa_project/from_ilifu/awigen_dataset/plink_file/awigen-qc"
## Keep only biallelic
${PLINK2} --bfile ${bfile} \
  --maf ${MINMAF} \
  --max-alleles 2 \
  --min-alleles 2 \
  --snps-only 'just-acgt' \
  --make-bed \
  --out ${out_prefix}/biallelic
# --maf ${MINMAF} 
# --max-alleles 2 
#  --snps-only ['just-acgt']

## LD prune
${PLINK2} --bfile ${out_prefix}/biallelic \
  --indep-pairwise 50 5 0.2 \
  --out ${out_prefix}/ld_pruned_biallelic

${PLINK2} --bfile ${out_prefix}/biallelic \
  --exclude ${out_prefix}/ld_pruned_biallelic.prune.out \
  --maf ${MINMAF} \
  --out ${out_prefix}/pruned \
  --recode ped \
  --geno 0.05

