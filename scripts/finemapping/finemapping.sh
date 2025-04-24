#!/bin/bash


set -e
FINEMAP="/home/yal084/tools/finemap_v1.4_x86_64/finemap_v1.4_x86_64"
config_n=50000
causal_n=10

infile_folder="/expanse/protected/gymreklab-dbgap/mount/H3Africa/DS_I_Africa_project/from_ilifu/results/finemapping"
all_infiles=$(cat ${infile_folder}/*/chr*infile*)

for file in ${all_infiles}; do 
  echo $(dirname ${file})/finemap_output.snp >> ${infile_folder}/finemapped_results_path.txt
  ${FINEMAP} --sss \
    --in-files ${file} \
    --log \
    --n-configs-top ${config_n} \
    --n-causal-snps ${causal_n}
done


