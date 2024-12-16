# TR imputation
## 1. Run LiftOver 
Run the `liftover.py` to update the POS from hg19 to hg38 and write lifted, unlifted SNPs into vcf files (2G memory). Use the `liftover.yaml` to set up the `pyliftover` package.  

The script will first remove SNPs with big AF difference between 1kg REF and target REF, then perform the liftover.

## 2. Run imputation and annotation
The TR imputation are based on this [TRTools](https://github.com/gymrek-lab/TRTools/blob/tr-gwas-tutorial/doc/VIGNETTE-GWAS-TUTORIAL.rst)

* Use the --window parameter to reduce memory useage (5 works for 10k sample with a 25G memory). A small window sometimes will have no marker SNPs inside, so use the default 40 and split samples to do the imputation.
```bash
Using --windowsize=5 
Running time for chromosome 21 
Using 
Haplotype phasing time:        3 hours 12 minutes 53 seconds
Imputation time:               1 hour 37 minutes 34 seconds
Total time:                    4 hours 52 minutes 31 seconds

Using --windowsize=40 (default one), takes about 6 hours to finish all
```
The extract and annotation part seems not memory intensive, 2G memory should be okay.
* The `'ID~"EnsTR"'` is not working, use `"ID=@file"` to include EnsTRs. `file` can be a text file contains all EnsTR IDs.
* Install `trtools` and use the `annotaTR` to add the annotations.  
