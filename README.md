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
## 1. Convert Plink to vcf format 
If SNPs in plink format, convert to vcf files using the scripts in `/scripts/plinkTovcf/plink2_conversion.sh`. This is pretty fast and can be done with 8G memory for 10k samples size. 

## 2. Imuptation
Before imputation, we need to check the quality of SNPs by comparing to the imputation reference panel, lift over the coordinate from hg19 to hg38, batch the samples to reduce the memory requirement by running:
```bash
# this step need about 2G memory
bash ./scripts/str_imputation/run_liftover.sh
```
Here are some filtering has been applied:
* remove SNPs where the AF difference is too big between 1kg (reference) and H3Africa cohort. 
* split the vcf files by samples, every 1,000 samples per batch.

After filtering and batching the vcf files, run the following scripts to impute STRs into SNPs: 
```bash 
# This need 25 GB memory for 1,000 samples
./scripts/str_imputation/beagle_imputation.sh
```

The last step is to annotate TRs and combine the samples by running:
```bash
## If only test TRs, run this scirpt to remove SNPs

## Run this scripts to include SNPs

```

