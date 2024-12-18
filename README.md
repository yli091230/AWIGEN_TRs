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
### 2.1 Prepare files for imputation
Before imputation, we need to check the quality of SNPs by comparing to the imputation reference panel, lift over the coordinate from hg19 to hg38, batch the samples to reduce the memory requirement by running:
```bash
# this step need about 2G memory
bash ./scripts/str_imputation/run_liftover.sh
```
Here are some filtering has been applied:
* remove SNPs where the AF difference is too big between 1kg (reference) and H3Africa cohort. 
* split the vcf files by samples, every 1,000 samples per batch.
### 2.2 Run imputation
After filtering and batching the vcf files, run the following scripts to impute STRs into SNPs: 
```bash 
# This need 25 GB memory for 1,000 samples
# To run in a interactive node
bash ./scripts/str_imputation/beagle_imputation.sh

# To submit jobs using SLURM
bash ./scripts/str_imputation/submit_beagle_jobs.sh
```
### 2.3 Merge and prepare files for GWAS
The last step is to combine the samples, annotate TRs and computing dosages using `annotaTR` from `TRtools` by running:
```bash
# chrom is a integer for chrom number
# str_only takes value "true" or "false". If true, will extract only TRs; if False, will include both TR and SNPs
bash ./scripts/str_imputation/extract_and_annotate_TRs.sh ${chrom} ${str_only}

# To submit jobs using SLURM
bash ./scripts/str_imputation/submit_annotatr.sh
```
It might be better to extract the TRs (by setting ${str_only}="true") and perform the GWAS for TR and SNPs seperately. The imputed variants may contains new SNPs. 


## 3. Running GWAS
GWAS on AWIGENE dataset have been reported and can be used as a reference [lipid traits](https://pmc.ncbi.nlm.nih.gov/articles/PMC9095599/), [blood pressure traits](https://www.nature.com/articles/s41467-023-44079-0#Sec10). 
