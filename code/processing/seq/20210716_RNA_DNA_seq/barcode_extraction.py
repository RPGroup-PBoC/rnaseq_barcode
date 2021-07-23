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
# Import this project's library
import rnaseq_barcode as rnaseq

# Find home directory for repo
repo = git.Repo("./", search_parent_directories=True)
homedir = repo.working_dir

# Date for sequencing run of the library to map
DATE = 20210716
# Description to be attached to folder names
DESCRIPTION = "_RNA_DNA_seq/"
# Define data directory
datadir = f"{homedir}/data/processed_sequencing/{DATE}{DESCRIPTION}"
# Read CSV file explaining what each file is
metadata = pd.read_csv(
    f"{homedir}/data/demux_sequencing/{DATE}{DESCRIPTION}/MANIFEST",
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
    # Extract barcode
    df_seq = df_seq.assign(op_bc=[x[0:20] for x in df_seq.sequence])

    # Define sequence next to GFP barcode
    ref_seq = "TTATTTGTACAGTT"

    # Initialize list to save index positions
    gfp_idx = list()
    gfp_bc = list()
    # Iterate over sequences and define the position of the barcode sequence
    for i, row in df_seq.iterrows():
        # Define regular expression to find sequence
        reg = regex.finditer(f"({ref_seq})" + "{e<=0}", row.sequence)
        # Extract the span for where the sequences was located
        match = [m.span() for m in reg]
        # Rules for lack of match
        if len(match) == 0:
            gfp_idx.append(np.nan)
            gfp_bc.append("NNNNN")
        # Rules for when there is a match
        else:
            gfp_idx.append(match[0][0] - 5)
            gfp_bc.append(row.sequence[match[0][0] - 5 : match[0][0]])

    df_seq = df_seq.assign(gfp_bc=gfp_bc, gfp_idx=gfp_idx)

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
