# Plink file notes

# conversion notes:
1. Need to provide reference allele of SNPs
    * can use plink2 based on this [link](https://groups.google.com/g/plink2-users/c/hCHiS9xXB5M)
    * The reference genome are downloaded from [1kg](http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/)
2. Check if the plink files are phased

# quality checks:
1. Use the 1GK AFR to check `/storage/resources/datasets/1000Genomes/phase3/`
2. Might need to remove "A/T", "C/G" ambigious SNPs
3. Using plink2 --missing to get the F_MISSING and filter on the F_MISSING. (This need ~8G memory for 10k samples)
