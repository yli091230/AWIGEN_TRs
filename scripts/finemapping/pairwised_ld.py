import pandas as pd
import sys
import numpy as np

snp_file = sys.argv[1]
tr_file = sys.argv[2]
variant_file = sys.argv[3]
out_dir = sys.argv[4]

'''
The dosage columns are named like ID_Allels e.g "19:123423:A:G_A"
need to remove the "_A" to keep the same as ID
'''


def format_name(v_name):
    split_name = v_name.split("_")
    if len(split_name) > 2:
        return "_".join(split_name[:-1])
    else:
        return split_name[0]


def format_df(dosage_df):
    variant_list = dosage_df.columns.tolist()[6:]
    dosage_df = dosage_df[["IID"]+variant_list].copy()
    dosage_df.columns = [format_name(i) for i in dosage_df.columns.tolist()]
    return dosage_df


# load dosage
snp_df = pd.read_csv(snp_file, sep="\t")
tr_df = pd.read_csv(tr_file, sep="\t")
variant_df = pd.read_csv(variant_file, header=None, names=["variants"])
variant_list = variant_df.variants.tolist()
# variant_list = [i.split("_")[0] for i in variant_list]

# combine dosage files and output the pearson r
snp_dosage = format_df(snp_df)
tr_dosage = format_df(tr_df)
merged_df = tr_dosage.merge(snp_dosage, on=["IID"]). set_index("IID")
merged_df = merged_df[variant_list].copy()
numpy_merged = merged_df.values
correlation_matrix = np.corrcoef(numpy_merged, rowvar=False)
np.savetxt(out_dir, correlation_matrix, delimiter=" ", fmt="%.5f")

# corr_df = merged_df.corr()
# # save to folder
# corr_df.to_csv(out_dir, sep=" ", index=False, header=False)
