# Plink file notes

# conversion notes:
1. Need to provide reference allele of SNPs
    * can use plink2 based on this [link](https://groups.google.com/g/plink2-users/c/hCHiS9xXB5M)
    * The reference genome are downloaded from [1kg](http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/)
2. Check if the plink files are phased

# quality checks:
Check the `../../notebooks/plinkVcf_quality_check/converted_vcf_quality_check.ipynb` for details.

1. Use the 1GK AFR to check `/storage/resources/datasets/1000Genomes/phase3/`
2. Might need to remove "A/T", "C/G" ambigious SNPs
* Ambigious SNPs looks not too bad, so they were kept in the imputation.
3. Using plink2 --missing to get the F_MISSING and filter on the F_MISSING. (This need ~8G memory for 10k samples) 
* F_MISSING tag is not used to filter SNPs for imputation.
4. Compare the Allele frequence:
* Filtering based on the Allele frequence (check the script for detailed method).
