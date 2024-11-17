# TR imputation
## LiftOver 
Run the `liftover.py` to update the POS from hg19 to hg38 and write lifted, unlifted SNPs into vcf files (2G memory).

## Imputation
The TR imputation are based on this [TRTools](https://github.com/gymrek-lab/TRTools/blob/tr-gwas-tutorial/doc/VIGNETTE-GWAS-TUTORIAL.rst)

Use the --window parameter to reduce memory useage (5 works for 10k sample with a 25G memory)
Running time for chromosome 21
Haplotype phasing time:        3 hours 12 minutes 53 seconds
Imputation time:               1 hour 37 minutes 34 seconds
Total time:                    4 hours 52 minutes 31 seconds

