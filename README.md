# AWI-GEN data
This the the dataset of `EGAD00010001996`. Plink files are downloaded from ilifu serve. Not sure how the file is generated but may worth to check this [method](https://github.com/h3abionet/h3agwas/tree/master/call2plink).  

The `EGAF00004691444` folder contains files sample informations, the `EGAF00004691445` contains description of each columns.

# File structure
Results and additional files are not included in this repository.

```bash
├── notebooks
├── README.md
├── results
└── scripts
    ├── plinkTovcf
    └── str_imputation
```
# Quick start
* Prepare files for imputation:
1. If SNPs in plink format, convert to vcf files and lift the coordinate to hg38 using the scripts in `plinkTovcf`.
2. Impute TRs (beagle is very memory intensive, for 10k samples, able to run imputation with 25G memory using `--window 5`). Extract TRs.
