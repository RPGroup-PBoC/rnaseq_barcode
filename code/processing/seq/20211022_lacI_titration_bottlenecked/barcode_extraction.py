#%%
import os
import glob
import itertools
import re
import regex
import numpy as np
import pandas as pd
import skbio
import collections
import git
import rnaseq_barcode as rnaseq

#%%

# Find home directory for repo
repo = git.Repo("./", search_parent_directories=True)
homedir = repo.working_dir

# Date for sequencing run of the library to map
DATE = 20211022
# Description to be attached to folder names
DESCRIPTION = 'lacI_titration_bottlenecked/'

# Find home directory for repo
repo = git.Repo("./", search_parent_directories=True)
homedir = repo.working_dir
datadir = homedir + f'/data/processed_sequencing/{DATE}{DESCRIPTION}'

# Find sequencing files
files = glob.glob(datadir + "*.fastq.gz")

# output directory
outdir = f"{homedir}/data/barcodes/{DATE}{DESCRIPTION}"

if not os.path.exists(outdir):
    os.makedirs(outdir)
    
# read reference GFP barcodes
df_gfp = pd.read_csv("./gfp_barcode.csv")


def find_close_bc(seq, barcodes):
    bc_lists = [list(x) for x in barcodes]
    for j, bc in enumerate(bc_lists):
        if np.sum([x != y for x, y in zip(list(seq), bc)]) == 1:
            return barcodes[j]
        
    return np.nan

def find_close_bc(seq, barcodes):
    bc_lists = [list(x) for x in barcodes]
    for j, bc in enumerate(bc_lists):
        if np.sum([x != y for x, y in zip(list(seq), bc)]) == 1:
            return barcodes[j]
        
    return np.nan

# %%
def analyze_file(file):
    print("Reading sequences from fastq file")
    df_seq = rnaseq.seq.read_sequences(file)

    # cDNA reads are treated differently
    RNA = False
    if "RNA" in file:
        RNA = True

    # Add index to dataframe
    df_seq["index"] = file.split("/")[-1]

    print("Mapping operator and GFP barcode")

    # Define sequence next to GFP barcode
    ref_seq_gfp = "TTATTTGTAC"
    ref_seq_op = "TAAATCCCACCCGATGCCTGCAGG"

    # Find GFP barcodes and filter for correct location
    gfp_bc_list = df_seq['sequence'].apply(rnaseq.seq.find_gfp_bc, args=(ref_seq_gfp, [50, 60], 1))
    df_seq['gfp_idx'], df_seq['gfp_bc'] = zip(*gfp_bc_list)
    op_bc_list = df_seq['sequence'].apply(rnaseq.seq.find_gfp_bc, args=(ref_seq_op, [0, 10], 1, 20, 'down'))
    df_seq['op_idx'], df_seq['op_bc'] = zip(*op_bc_list)

    if RNA:
        df_seq['UMI'] = [x[0:10] for x in df_seq.sequence]

    
    print("Filtering GFP barcodes")
    df_filt = df_seq.dropna(subset=['gfp_bc'])

    # Define if barcode is present
    bc_bool = [x not in df_gfp.gfp_barcode_revcomp.values for x in df_filt.gfp_bc.values]

    df_filt.loc[bc_bool,'gfp_bc'] = df_filt.loc[bc_bool,'gfp_bc'].apply(find_close_bc, args=(df_gfp.gfp_barcode_revcomp.values,))
    df_filt = df_filt.dropna(subset=['gfp_bc'])

    # Keep only desired columns
    df_filt = df_filt[["sequence", "index", "op_bc", "gfp_bc"]]
    print("Saving into memory")

    file_name = file.split("/")[-1].replace('.fastq.gz', '.csv')
    df_filt.to_csv(f"{outdir}{file_name}", index=False)

    
# Parallelize analysis    
a_pool = multiprocessing.Pool(cores)
result = a_pool.map(analyze_file, files)
