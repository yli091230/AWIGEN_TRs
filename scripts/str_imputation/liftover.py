#!/usr/bin/env python3

import gzip
import pandas as pd
from io import StringIO
from pyliftover import LiftOver
import sys
import logging
import numpy as np

chain_file="/expanse/protected/gymreklab-dbgap/mount/H3Africa/DS_I_Africa_project/from_ilifu/additional_files/liftover_file/hg19ToHg38.over.chain.gz"
lo = LiftOver(chain_file)


# Function to read VCF.gz file into a pandas DataFrame
def read_vcf_gz(file):
    """
    read the input vcf.gz file
    output the header and SNPs info in pandas dataframe format
    """
    with gzip.open(file, 'rt') as f:
        # Extract header and body lines
        header_lines = []
        body_lines = []
        for line in f:
            if line.startswith('##'):
                header_lines.append(line)
            elif line.startswith('#CHROM'):
                header_lines.append(line)
                body_lines.append(line)
            else:
                body_lines.append(line)
        # Create DataFrame from body lines
        df = pd.read_csv(StringIO(''.join(body_lines)), sep='\t')
    return header_lines, df


def save_vcf(df, header_lines, output_file):
    """
    write the head_lines and df to vcf.gz files
    """
    with open(output_file, 'wt') as f:
        # Write header lines
        for line in header_lines:
            f.write(line)
        # Write DataFrame to VCF format
        df.to_csv(f, sep='\t', index=False, header=False)


def hg19Tohg38_batch(chrom, pos):
    """
    Batch liftover for hg19 to hg38.
    Returns an array of converted positions, or None if lift failed.
    """
    result = []
    for c, p in zip(chrom, pos):
        try:
            lifted = lo.convert_coordinate(c, p)
            if lifted:
                result.append(lifted[0][1])
            else:
                result.append(None)
        except (IndexError, TypeError):
            result.append(None)
    return np.array(result)


def hg19Tohg38(hg19_chr, hg19_pos):
    """
    liftover hg19 to hg38
    return None if lift failed
    """
    try:
        return lo.convert_coordinate(hg19_chr, hg19_pos)[0][1]
    except (IndexError, TypeError):
        return None


#def read_vcf_gz_chunked(input_vcf, chunksize=10000):
#    """
#    Read VCF file in chunks
#    """
#    header_lines = []
#    with open(input_vcf, 'rt') as f:
#        for line in f:
#            if line.startswith('#'):
#                header_lines.append(line)
#            else:
#                break
#    return header_lines, pd.read_csv(input_vcf, comment='##', delim_whitespace=True, chunksize=chunksize, compression="gzip")


def liftCoordinates_chunked(input_vcf, logger=None):
    """
    input_vcf is on hg19, check if the chrom starts with "chr" if not add back.
    return two dataframes: the lifted and unlifted ones
    """
    header_lines, chunk = read_vcf_gz(input_vcf)
    df_lifted_list = []
    df_unlifted_list = []

#    for chunk in chunks:
    if logger is not None:
        logger.info("Processing chunk...")
    # Add "chr" prefix if missing
    if "chr" not in str(chunk.loc[0]["#CHROM"]):
        chunk["#CHROM"] = "chr" + chunk["#CHROM"].astype(str)
    # Perform batch liftover
    chunk["NEW_POS"] = hg19Tohg38_batch(chunk["#CHROM"].values, chunk["POS"].values)
    # Separate lifted and unlifted rows
    df_lifted_chunk = chunk[~chunk["NEW_POS"].isnull()].copy()
    df_lifted_chunk["POS"] = df_lifted_chunk["NEW_POS"].astype(int)
    df_lifted_chunk = df_lifted_chunk.drop(columns=["NEW_POS"])
    df_unlifted_chunk = chunk[chunk["NEW_POS"].isnull()].copy()
    df_unlifted_chunk = df_unlifted_chunk.drop(columns=["NEW_POS"])
    df_lifted_list.append(df_lifted_chunk)
    df_unlifted_list.append(df_unlifted_chunk)

    df_lifted = pd.concat(df_lifted_list, ignore_index=True)
    df_unlifted = pd.concat(df_unlifted_list, ignore_index=True)
    return df_lifted, df_unlifted, header_lines


def liftCoordinates(input_vcf, logger=None):
    """
    input_vcf is on hg19, check if the chrom starts with "chr" if not add back.
    return two datafram: the lifted and unlifted ones
    """
    header_lines, df = read_vcf_gz(input_vcf)
    if logger is not None:
        logger.info("Finish loading, start liftOver ...")
    if "chr" not in str(df.loc[0]["#CHROM"]):
        df["#CHROM"] = "chr"+df["#CHROM"].astype(str)
    df["New_POS"] = df[["#CHROM", "POS"]].apply(lambda x: hg19Tohg38(x["#CHROM"], x["POS"]), axis=1)
    df_lifted = df[~df["New_POS"].isnull()].copy()
    df_lifted["POS"] = df_lifted["New_POS"].astype(int)
    df_lifted.drop(columns="New_POS", inplace=True)
    df_unlifted = df[df["New_POS"].isnull()].copy()
    df_unlifted.drop(columns="New_POS", inplace=True)

    return df_lifted, df_unlifted, header_lines


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

    df_hg38, df_failed, headers = liftCoordinates(input_file, logger=logger)
    logger.info("Finish liftover, and start to save files ...")
    save_vcf(df_hg38, headers, f"{output_dir}/{vcf_name}_hg38.vcf")
    save_vcf(df_failed, headers, f"{output_dir}/{vcf_name}_unlifted.vcf")
    logger.info(f"Files saved. There are {len(df_hg38):,} lifted SNPs, {len(df_failed):,} SNPs failed to lift.")
