#!/usr/bin/env python3

import gzip
from pyliftover import LiftOver
import sys
import logging

# can specify chain file or use liftOver("hg19", "hg38") to download online automatically
chain_file="/expanse/protected/gymreklab-dbgap/mount/H3Africa/DS_I_Africa_project/from_ilifu/additional_files/liftover_file/hg19ToHg38.over.chain.gz"
lo = LiftOver(chain_file)


def update_vcf_gz(input_file, output_lifted, output_failed):
    """
    read the input vcf.gz file
    liftover POS line by line
    """
    lifted_num = 0
    failed_num = 0
    with gzip.open(input_file, 'rt') as infile, open(output_lifted, 'w') as outfile_lifted, open(output_failed, "w") as outfile_failed:
        for line in infile:
            if line.startswith('#'):
                outfile_lifted.write(line)
                outfile_failed.write(line)
            else:
                vcf_values = line.split("\t")
                curr_chrom = vcf_values[0]
                curr_pos = int(vcf_values[1])
                curr_chrom = curr_chrom if curr_chrom.startswith("chr") else f"chr{curr_chrom}"
                lifted_pos = hg19Tohg38(curr_chrom, curr_pos)
                if lifted_pos is not None:
                    vcf_values[1] = str(lifted_pos)
                    new_line = "\t".join(vcf_values)
                    outfile_lifted.write(new_line)
                    lifted_num += 1
                else:
                    outfile_failed.write(line)
                    failed_num += 1
    return lifted_num, failed_num


def hg19Tohg38(hg19_chr, hg19_pos):
    """
    liftover hg19 to hg38
    return None if lift failed
    """
    try:
        return lo.convert_coordinate(hg19_chr, hg19_pos)[0][1]
    except (IndexError, TypeError):
        return None


if __name__ == "__main__":

    input_file = sys.argv[1]
    output_dir = sys.argv[2]
    vcf_name = input_file.split("/")[-1].split(".")[0]
    log_file = f"{output_dir}/{vcf_name}_liftOver.log"
    logging.basicConfig(filename=log_file,
                        format='%(asctime)s %(message)s',
                        filemode="w")
    logger = logging.getLogger()
    logger.setLevel(logging.INFO)
    logger.info(f"start loading {vcf_name}.vcf.gz ...")
    lifted_num, failed_num = update_vcf_gz(input_file, f"{output_dir}/{vcf_name}_hg38.vcf", f"{output_dir}/{vcf_name}_unlifted.vcf")
    logger.info(f"Finish liftover: {lifted_num:,} lifted to hg38, {failed_num:,} SNPs failed.")
