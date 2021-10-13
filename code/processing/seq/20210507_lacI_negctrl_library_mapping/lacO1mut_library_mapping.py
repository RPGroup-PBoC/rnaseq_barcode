#%%
import os
import glob
import itertools
import re
import numpy as np
import pandas as pd
import collections
import skbio
import git

#%%
# Find project parental directory
repo = git.Repo("./", search_parent_directories=True)
homedir = repo.working_dir

# Define data directory
datadir = f"{homedir}/data/processed_sequencing/20210507_lacI_negctrl_library_mapping/"

# Define output dir
outputdir = f"{homedir}/data/barcodes/20210507_lacI_negctrl_library_mapping/"

# List fastq.gz file
fastq_file = glob.glob(f"{datadir}*LacO1_mutants*.fastq.gz")[0]

#%%
# O1 binding site sequence
O1_seq = "aattgtgagcggataacaatt"

# Forward primer
fwd_prim = skbio.DNA("GCTTATTCGTGCCGTGTTAT").reverse_complement()
# Reverse primer
rev_prim = skbio.DNA("GGGCACAGCAATCAAAAGTA").reverse_complement()

# Define RNAP binding site
rnap = str(
    skbio.DNA("TTTACACTTTATGCTTCCGGCTCGTATAATGTGTGG").reverse_complement()
)

# Define clone binding site
clone = str(skbio.DNA("gctagcCAATGCGGgagctc".upper()).reverse_complement())

#%%
print("Reading lacO1mut sequences into memory")

# Use skbio to have a generator to iterate over fastq
seqs = skbio.io.read(
    f"{fastq_file}", format="fastq", verify="false", variant="illumina1.8"
)

# Initialize list to save sequences
seq_list = list()

# Initialize counter
counter = 0

# Define number of samples
n_samples = 10000

# Iterate over sequences
for seq in seqs: #itertools.islice(seqs, n_samples):
    if counter % 10000 == 0:
        print(f"reading seq #{counter}")
    # Extract sequence information
    seq_id = seq.metadata["id"]
    sequence = str(skbio.DNA(sequence=seq, validate=False))
    # Append to list
    seq_list.append([seq_id, sequence])

    counter += 1

# Initialize dataframe to save sequences
names = ["id", "sequence"]
df_seq = pd.DataFrame.from_records(seq_list, columns=names)

# Add index and sequence length to dataframe
df_seq["seq_len"] = df_seq["sequence"].apply(len)

print("Done reading sequences...")

#%%

print("Mapping forward and reverse primers")
# Initialize array to save primer start position
prim_pos = np.zeros([len(df_seq), 2], dtype=int)

# Loop through sequences
for i, seq in df_seq.iterrows():
    # Search forward primer
    fwd_pos = re.search(str(fwd_prim), seq["sequence"])
    # Save position
    if bool(fwd_pos):
        prim_pos[i, 0] = fwd_pos.span()[0]

    # Search reverse primer
    rev_pos = re.search(str(rev_prim), seq["sequence"])
    # Save position
    if bool(rev_pos):
        prim_pos[i, 1] = rev_pos.span()[0]

# Assing columns with information
df_seq = df_seq.assign(
    fwd_prim=prim_pos[:, 0],
    rev_prim=prim_pos[:, 1],
    prim_dist=np.abs(prim_pos[:, 0] - prim_pos[:, 1]),
)
print("Done mapping primers...")

#%%
# Filtering sequences
print("Filtering sequences by:")

print("1. Filtering by separation between primers")

# Save original dataframe length
len_df = len(df_seq)

# Filter by length
df_filt = df_seq[df_seq["prim_dist"] == 82]

# Save rejected sequences
df_reject = df_seq[df_seq["prim_dist"] != 82][["id", "sequence"]]
df_reject = df_reject.assign(reject_by="primers_distance")

# Reset index
df_filt.reset_index(inplace=True, drop=True)
df_reject.reset_index(inplace=True, drop=True)

# Print percentage of sequences removed
print(
    f"""
Cumulative percentage of original sequences removed: 
{100 - np.round(len(df_filt) / len_df * 100, 2)}%
"""
)

#%%
print("2. Filtering by RNAP binding site sequence and position")

# Initialize array to save RNAP position
rnap_pos = np.zeros(len(df_filt))
# Loop through sequences
for i, seq in df_filt.iterrows():
    # Search RNAP sequence
    rnap_re = re.search(str(rnap), seq["sequence"])
    # Save position
    if bool(rnap_re):
        rnap_pos[i] = int(rnap_re.span()[0])

# Add column to dataframe
df_filt = df_filt.assign(rnap=rnap_pos)

# Compute RNAP distance
rnap_dist = (
    (df_filt["rnap"] - df_filt["rev_prim"])
    - len(rev_prim)
    - len(O1_seq)
)

# Copy filtered step
df = df_filt

# Select sequences
df_filt = df[rnap_dist == 0]
# Store rejected sequences
df_r = df[rnap_dist != 0][["id", "sequence"]]
df_r = df_r.assign(reject_by="rnap_pos")
df_reject = df_reject.append(df_r, ignore_index=True)

# Reset index
df_filt.reset_index(inplace=True, drop=True)
df_reject.reset_index(inplace=True, drop=True)

# Print percentage of sequences removed
print(
    f"""
Cumulative percentage of original sequences removed: 
{np.round(100 - len(df_filt) / len_df * 100, 2)}%
"""
)
#%%
print("3. Filtering by cloning site sequence and position")

# Initialize array to save clone position
clone_pos = np.zeros(len(df_filt))
# Loop through sequences
for i, seq in df_filt.iterrows():
    # Search clone sequence
    clone_re = re.search(str(clone), seq["sequence"])
    # Save position
    if bool(clone_re):
        clone_pos[i] = clone_re.span()[0]

# Add column to dataframe
df_filt = df_filt.assign(clone=clone_pos)

# Compute clone distance
clone_dist = (df_filt["clone"] - df_filt["rev_prim"]) + len(clone)

# Copy filtered step
df = df_filt

# Select sequences
df_filt = df[(clone_dist == 0) & (clone_pos == 20)]
# Store rejected sequences
df_r = df[(clone_dist != 0) | (clone_pos != 20)][["id", "sequence"]]
df_r = df_r.assign(reject_by="clone_pos")
df_reject = df_reject.append(df_r, ignore_index=True)

# Reset index
df_filt.reset_index(inplace=True, drop=True)
df_reject.reset_index(inplace=True, drop=True)

# Print percentage of sequences removed
print(
    f"""
Cumulative percentage of original sequences removed: 
{np.round(100 - len(df_filt) / len_df * 100, 2)}%
"""
)

#%%

# Initialize dataframe to save outcome
# names = ["sequence", "barcode", "counts"]
# df_counts = pd.DataFrame(columns=names)

# Count unique sequences
df_op = df_filt["sequence"].value_counts().reset_index(name="counts")
df_op.columns = ["sequence", "counts"]


# Extract barcodes and operator
df_op = df_op.assign(
    barcode=df_op["sequence"].apply(lambda x: x[0:20]),
    operator=df_op["sequence"].apply(lambda x: x[81-len(O1_seq):81]),
)

#%%
# Write file to memory
print("writing barcode list into memory")
df_op.to_csv(f"{outputdir}lacO1mut_barcodes_counts.csv", index=False)

df_reject.to_csv(f"{outputdir}lacO1mut_rejected_seq.csv", index=False)

print("Done! Barcodes filtered and quantified")
# %%
