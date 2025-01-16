#!/bin/bash

: <<'comment'
  Need:
    PLINK2, python3
  This script take regression results from SNP and TRs, then
  1. Find variants with p value < 5e-8
  2. Select all nominal significant variants (p <0.05) around 100kb of that GWAS signal.
  3. Using those selected variants for finemapping

  The ld is calculated in python script
comment

chrom=$1
set -xe

gt_folder="/expanse/protected/gymreklab-dbgap/mount/H3Africa/DS_I_Africa_project/from_ilifu/results/imputed_TRs"
reg_folder="/expanse/protected/gymreklab-dbgap/mount/H3Africa/DS_I_Africa_project/from_ilifu/results/GWAS/associations/friedewald_ldl_c"

snp_reg=${reg_folder}/merged_chr${chrom}_SNPs_gwas.friedewald_ldl_c_c_qc.glm.linear
tr_reg=${reg_folder}/merged_chr${chrom}_TRs_gwas.friedewald_ldl_c_c_qc.glm.linear
snp_pfile=${gt_folder}/chr${chrom}/chr${chrom}_imputed_filtered_SNPs
tr_pfile=${gt_folder}/chr${chrom}/chr${chrom}_imputed_and_annotated_TRs


out_dir="/expanse/protected/gymreklab-dbgap/mount/H3Africa/DS_I_Africa_project/from_ilifu/results/finemapping/chr${chrom}"

PLINK2="/expanse/protected/gymreklab-dbgap/mount/H3Africa/DS_I_Africa_project/from_ilifu/scripts/plinkTovcf/plink_2.0/plink2"
CALCULATE_LD="/expanse/protected/gymreklab-dbgap/mount/H3Africa/DS_I_Africa_project/from_ilifu/scripts/finemapping/pairwised_ld.py"
#FINEMAP="/home/yal084/tools/finemap_v1.4_x86_64/finemap_v1.4_x86_64"
### the out_dir should contains the chromsome already
finemapping_dir=${out_dir}/finemapping_loci
mkdir -p ${finemapping_dir}

# default 1M bp
flnk_range=1000000

filtering_snps() {
  file_name=$1
  p_threshold=$2
  out_file=$3
  first=$4
  if [ "$first" -eq 0 ]; then
    # select variatns pass the p_threshold and have no error code from the regression
    awk -F"\t" -v p_val=${p_threshold} '$(NF-1)+0<p_val && $NF=="."' ${file_name} > ${out_file}
  else 
    awk -F"\t" -v p_val=${p_threshold} '$(NF-1)+0<p_val && $NF=="."' ${file_name} >> ${out_file}
  fi

}
echo "Extract sig (nominal or 5e-8) gwas loci"
## 1.extract all gwas hits
gwas_hit=${out_dir}/gwas_hit.tsv
gwas_hit_sorted=${out_dir}/gwas_hit_sorted.tsv

filtering_snps "${snp_reg}" "5e-8" ${gwas_hit} 0
filtering_snps "${tr_reg}" "5e-8" ${gwas_hit} 1

sort -k2,2n -u ${gwas_hit} > ${gwas_hit_sorted} 
rm ${gwas_hit}

# extract all nominal significant
nominal_sig=${out_dir}/nominal_sig.tsv
nominal_sig_sorted=${out_dir}/nominal_sig_sorted.tsv

filtering_snps "${snp_reg}" "0.05" ${nominal_sig} 0
filtering_snps "${tr_reg}" "0.05" ${nominal_sig} 1
# sort variants by p-values
sort -k2,2n -u ${nominal_sig} > ${nominal_sig_sorted}
rm ${nominal_sig}
echo "Finished extraction, search for leading GWAS loci"
## 2.get list of locis for fine-mapping
gwas_loci=${out_dir}/gwas_loci.tsv
awk -F"\t" -v window=2000000 -v subwindow=1000000 '
    BEGIN {
        lower_boundary = 0;
        selected_boundary = 0;
        last_boundary = 0;
        out_line = "";
    }
    {
        if (NR == 1) {
            print $0;
            next;
        }
        # Initialize the first entry
        if (NR == 2) {
            lower_boundary = $2;
            selected_boundary = $2;
            old_boundary = $2;
            out_line = $0;
            next;
        }

        # Check if within subwindow
        if ($2 <= lower_boundary + subwindow) {
            old_boundary = selected_boundary
            selected_boundary = $2;
            out_line = $0;
        }
        else if ($2 <= lower_boundary + window) {
            # Within the 2MB window but beyond the subwindow
            if ($2 < old_boundary + subwindow) {
                selected_boundary = $2;
            }
            else {
                print out_line
                lower_boundary = $2;
                last_boundary = old_boundary
                out_line = $0;
            }
        }
        else {
            # Beyond 2MB window, finalize the region
            print out_line;
            last_boundary = old_boundary
            lower_boundary = $2;
            selected_boundary = $2;
            out_line = $0;
        }
    }
    END {
        # Print last region if applicable
        if ($2 > last_boundary + subwindow) {
            print out_line;
        }
    }' ${gwas_hit_sorted} > ${gwas_loci}
echo "Finish loci identification, start preparation FINEMAP files..."
finemapping_infiles=${out_dir}/finemappling_infile_list.tsv
## 3.prepare finemapping files for each locus
while IFS= read -r line; do
  # skipping the header line
  if [[ ${line} == \#* ]]; then
    continue
  fi
  IFS=$'\t' read -r -a v_array <<< ${line}
  variant_pos=${v_array[1]}
  variant_id=${v_array[2]}
  sample_size=${v_array[10]}
  lower_boundary=$(( ${variant_pos} - ${flnk_range} ))
  upper_boundary=$(( ${variant_pos} + ${flnk_range} ))
  # create a folder for current finemapping region
  curr_folder=${finemapping_dir}/chr${chrom}:${lower_boundary}-${upper_boundary}
  echo "    Processing $(basename ${curr_folder})" 
  mkdir -p ${curr_folder}
  
  # output the master file
  echo  "z;ld;snp;config;cred;log;n_samples" > ${curr_folder}/finemap_input.master
  echo ${curr_folder}/{finemap_input.z,finemap_input.ld,finemap_output.snp,finemap_output.config,finemap_output.cred,finemap_output.log,${sample_size}} | sed 's/ /;/g' >> ${curr_folder}/finemap_input.master
  echo "    output z file...)"
  # output the z-file
  awk -F"\t" -OF" " -v low_b=${lower_boundary} -v up_b=${upper_boundary} '$2+0<up_b && $2+0>low_b {print $3,$1,$2,"nan","nan","nan",$12,$13}' ${nominal_sig_sorted} > ${curr_folder}/finemap_input_no_sorted.z
  sort -k3,3n ${curr_folder}/finemap_input_no_sorted.z > ${curr_folder}/finemap_input.z
  
  ## prepare the ld file
  # get the variants ID 
  awk -F" " 'NR>1 {print $1}' ${curr_folder}/finemap_input.z > ${curr_folder}/variants_for_finemapping.txt
  ${PLINK2} --pfile ${tr_pfile} \
    --extract ${curr_folder}/variants_for_finemapping.txt \
    --export A \
    --out ${curr_folder}/tr_dosage

  ${PLINK2} --pfile ${snp_pfile} \
    --extract ${curr_folder}/variants_for_finemapping.txt \
    --export A \
    --out ${curr_folder}/snp_dosage
   echo "    calculating LD...)" 
  
  python3 ${CALCULATE_LD} ${curr_folder}/snp_dosage.raw \
    ${curr_folder}/tr_dosage.raw \
    ${curr_folder}/variants_for_finemapping.txt \
    ${curr_folder}/finemap_input.ld
  # create a file to store all finemapping master files
  echo ${curr_folder}/finemap_input.master >> ${finemapping_infiles}
done < ${gwas_loci}


# ## extract corresponding SNPs and STRs
#
#
# ## combine SNPs and STRs pgen file
# ${PLINK2} \
# 	--pmerge-list-dir ${curr_signals} \
#   --out
#
# ## calculate LD files
# ${PLINK2} \
# 	--bfile /data/module1/Ref_bfile/g1000_eur \
# 	--keep-allele-order \
# 	--ld-window-r2 0 \
# 	--r-phased \
# 	--extract locus1_snps.txt \
# 	--out locus1_1000G
