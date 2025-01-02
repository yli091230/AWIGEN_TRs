# AWI-GEN data
This the the dataset of `EGAD00010001996`. Plink files are downloaded from ilifu serve. Not sure how the file is generated but may worth to check this [method](https://github.com/h3abionet/h3agwas/tree/master/call2plink).  

The `EGAF00004691444` folder contains files sample informations, the `EGAF00004691445` contains description of each columns.

# File structure
Results and additional files are not included in this repository.

```bash
├── additional_files
│   ├── AFR_AF_from_1kg
│   ├── dbSNP
│   ├── genetic_maps
│   ├── h3a_A3
│   ├── liftover_file
│   └── reference_genome
├── awigen_dataset
│   ├── phenotypes
│   ├── plink_file
│   └── WGS_TR_call_from_Ibra
├── notebooks
│   ├── gwas_results
│   ├── Phenotypes_check
│   └── plinkVcf_quality_check
├── other_gwas_results
│   └── margoliash-et-al-2023
├── results
│   ├── converted_lifted_vcf
│   ├── converted_vcf
│   ├── GWAS
│   └── imputed_TRs
└── scripts
    ├── GWAS
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

### Use PLINK2 for GWAS association
There are different tools avaiable for GWAS association test. Here we use PLINK2 to run GWAS.

### 3.1.1 Covariates for GWAS  
GWAS is performed seperately for SNPs and TRs. The most common used covariates contains: age, sex, PCs. Dependents on phenotypes, other covariates or adjustement on pehnotypes need to be added. Make sure to check literature on how those GWAS are performed for the specific phenotypes. In this example, LDL phenotype is used for GWAS and cholesterol treatment is included in the covariates.  
### PCs
To calcualte PCs, the input SNPs were first fitered to remove SNP with MAF < 0.05, and keeps only biallelic SNPs, then LD pruned. The filtered SNPs were then used for PCA calculation.
```bash
# To filtering SNPs
bash ./scripts/GWAS/prepration/filtering_SNPs.sh
# Use smartpca to calcualte PCA
bash ./scripts/GWAS/prepration/calculate_pc.sh 
```
### 3.1.2 QC on variants
Can use the filtered SNPs from last step for GWAS test.

### 3.1.3 Sample relateless check
It has been shown that some samples from the AWI-GEN dataset are closely related. PLINK2 doesn't support LMM model, but it do provide method to remove related samples. Use the --make-king to generate kinship and set the cutoff to 0.177

Tried to calculate the kinship use both all raw downloaded QCed SNPs or filtered one for PC. There are no difference. 
### 3.1.4 Run GWAS with PLINK2
* SNPs showing missingness greater than 0.05, MAF less than 0.01, and HardyWeinberg equilibrium (HWE) P-value less than 0.0001 are removed (the SNPs get from EGA should be QCed).
* Duplicates, sexual chromosomes, mitochondrial SNPs and SNPs failed to match the reference alleles are also removed (should be removed in the QC). 
* `"--mac 20"` is a reasonable filter  to apply before --glm



