# Finemapping


There are lots of different finemaping methods available. We use [FINEMAP](http://www.christianbenner.com) as an example.

## 1. Selection of finemapping range
Finemapping range can be selected using LD or fixed range of variant within the leading varaints. Here, we choose variants within 1Mb upstrean/downstream of an GWAS hit (p value < 5e-8, for details check the `./define_locus.sh`).

## 2. Input variants for finemapping
We use all variants with nominal p<0.05 for fine-mapping.

## 3. Number of causal variants
Number of causal variants need to try. Based on experience, can use the number where the PIP for max number of causal variants is 0.
