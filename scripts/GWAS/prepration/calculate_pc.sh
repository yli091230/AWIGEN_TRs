#!/usr/bin/env bash

export LD_LIBRARY_PATH=/home/yal084/.conda/envs/eign/lib:$LD_LIBRARY_PATH
export PATH=/expanse/protected/gymreklab-dbgap/mount/H3Africa/DS_I_Africa_project/from_ilifu/scripts/GWAS/prepration/EIG-7.2.1/bin:/expanse/protected/gymreklab-dbgap/mount/H3Africa/DS_I_Africa_project/from_ilifu/scripts/GWAS/prepration/EIG-7.2.1/src/eigensrc/:$PATH
CONVERTF="/expanse/protected/gymreklab-dbgap/mount/H3Africa/DS_I_Africa_project/from_ilifu/scripts/GWAS/prepration/EIG-7.2.1/src/convertf"
SMARTPCA="/expanse/protected/gymreklab-dbgap/mount/H3Africa/DS_I_Africa_project/from_ilifu/scripts/GWAS/prepration/EIG-7.2.1/bin/smartpca.perl"

PCAINPUT="/expanse/protected/gymreklab-dbgap/mount/H3Africa/DS_I_Africa_project/from_ilifu/results/GWAS/PCs/filtering_SNPs"
OUTPATH="/expanse/protected/gymreklab-dbgap/mount/H3Africa/DS_I_Africa_project/from_ilifu/results/GWAS/PCs/calculate_pc"

parfile="${OUTPATH}/convertf_parfile.txt"
#echo "genotypename: " ${PCAINPUT}/pruned.ped > ${parfile}
#echo "snpname: " ${PCAINPUT}/pruned.map >> ${parfile}
#echo "indivname: " ${PCAINPUT}/pruned.ped >> ${parfile}
#echo "outputformat: EIGENSTRAT" >> ${parfile}
#echo "genooutfilename: " ${OUTPATH}/eigenstratgeno >> ${parfile}
#echo "snpoutfilename: " ${OUTPATH}/snp >> ${parfile}
#echo "indoutfilename: " ${OUTPATH}/ind >> ${parfile}
#echo "familynames: NO" >> ${parfile}
#${CONVERTF} -p ${parfile}

# https://github.com/DReichLab/EIG/blob/master/EIGENSTRAT/README#L32-L65
${SMARTPCA} -i ${OUTPATH}/eigenstratgeno \
  -a ${OUTPATH}/snp \
  -b ${OUTPATH}/ind_country \
  -k 10 \
  -o ${OUTPATH}/pca \
  -e ${OUTPATH}/evals \
  -p ${OUTPATH}/plot \
  -l ${OUTPATH}/log

