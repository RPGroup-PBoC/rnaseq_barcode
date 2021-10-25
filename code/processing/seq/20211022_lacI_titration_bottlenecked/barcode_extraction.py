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

#%%

# Find home directory for repo
repo = git.Repo("./", search_parent_directories=True)
homedir = repo.working_dir

# Date for sequencing run of the library to map
DATE = 20211022
# Description to be attached to folder names
DESCRIPTION = 'lacI_titration_bottlenecked/'
# Define data directory
datadir = f"{homedir}/data/processed_sequencing/{DATE}_{DESCRIPTION}"
# Read CSV file explaining what each file is
metadata = pd.read_csv(
    f"{homedir}/data/demux_sequencing/{DATE}_{DESCRIPTION}/MANIFEST",
    comment="#"
)
metadata = metadata[metadata.direction == "forward"].reset_index()
metadata = metadata.rename(columns={"sample-id": "id"})

# output directory
outdir = f"{homedir}/data/barcodes/{DATE}{DESCRIPTION}/"

# read reference GFP barcodes
df_gfp = pd.read_csv("./gfp_barcode.csv")

#%%

for i, f in enumerate(metadata.filename):
    print(f)

    print("loading fastq into memory")
    # Use skbio to have a generator to iterate over fastq
    seqs = skbio.io.read(
        f"{datadir}{f}", format="fastq", verify="false", variant="illumina1.8"
    )

    # Initialize list to save sequence objects
    seq_list = list()
    # Initialize counter
    counter = 0

    # Iterate over sequences
    for seq in seqs:
        if counter % 100000 == 0:
            print(f"read # {counter}")
        # Extract sequence information
        seq_id = seq.metadata["id"]
        # Extract sequence
        sequence = str(skbio.DNA(sequence=seq, validate=False))
        # Append to list
        seq_list.append([seq_id, sequence])
        counter += 1

    # Initialize dataframe to save sequences
    names = ["id", "sequence"]
    df_seq = pd.DataFrame.from_records(seq_list, columns=names)

    # Add index and sequence length to dataframe
    df_seq["index"] = f
    df_seq["seq_len"] = df_seq.sequence.apply(len)

    print("Mapping operator and GFP barcode")

    # Define sequence next to GFP barcode
    ref_seq = "TAAATCCCACCCGATG"

    # Initialize list to save index positions
    gfp_idx = list()
    gfp_bc = list()
    op_bc = list()

    RNA = False
    if "RNA" in metadata.filename:
        RNA = True
        UMI_list = list()

    # Iterate over sequences and define the position of the barcode sequence
    for i, row in df_seq.iterrows():
        # Define regular expression to find sequence 
        reg = regex.finditer(f"({ref_seq})" + "{e<=0}", row.sequence)
        # Extract the span for where the sequences was located
        match = [m.span() for m in reg]
        # Rules for lack of match
        if len(match) == 0:
            gfp_idx.append(np.nan)
            gfp_bc.append("NNNN")
            op_bc.append(20*"N")
        # Rules for when there is a match
        else:
            gfp_idx.append(match[0][1] + 34)
            gfp_bc.append(row.sequence[match[0][1] +34: match[0][1]+39])
            op_bc.append(row.sequence[match[0][1] +8: match[0][1]+29])
            if RNA:
                UMI_list.append(row.sequence[0:10])

    if RNA:
        df_seq = df_seq.assign(gfp_bc=gfp_bc, gfp_idx=gfp_idx, op_bc=op_bc, UMI=UMI_list)
    else:
        df_seq = df_seq.assign(gfp_bc=gfp_bc, gfp_idx=gfp_idx, op_bc=op_bc)

    print("Applying barcode length and position filter")
    # Filter by position of GFP barcode
    df_seq = df_seq[df_seq.gfp_idx == 26]
    # Define if barcode is present
    bc_bool = [
        x in df_gfp.gfp_barcode_revcomp.values for x in df_seq.gfp_bc.values
    ]
    df_seq = df_seq[bc_bool]

    # Keep only desired columns
    df_seq = df_seq[["sequence", "index", "op_bc", "gfp_bc"]]
    print("Saving into memory")
    # Generate output directory if it doesn't exist
    if not os.path.exists(outdir):
        os.makedirs(outdir)

    df_seq.to_csv(f"{outdir}{f.replace('.fastq.gz', '.csv')}", index=False)
