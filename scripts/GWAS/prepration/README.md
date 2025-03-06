# GWAS file prepration

1. SNP QC and PCA calculation was performed as described in this [paper](https://www.nature.com/articles/s41467-022-30098-w).

* PCA was calculated using [EIGENSOF](https://hsph.harvard.edu/research/price-lab/software/)
* To plot by country of samples, the third columns of ind were replaced by country information.

2. Phenotypes and covariates prepration 

Here is an example of [notebook](../../../notebooks/Phenotypes_check/Phenotypes_check.ipynb) on prepreation

For covariates, sex and age information were extracted and remove samples with missing info. Joint with PCs to generate the covariates table. Other trait specific covariates may reuired, check literature for details.

Phenotypes are extracted from the `EGA_dataset_v0_1.csv` and perform the following filtering:
* Remove invalid values like `Not Applicable`, `Missing`, `out-of-range`, etc values based on the description of those measurment.
* Some phenotypes may need to be corrected, make sure to keep it consistent with literature.
* Save the phenotype in a compatiable PLINK2 format, which contains two columns: use the `study_id` as `#IID` and another column named as the trait.
