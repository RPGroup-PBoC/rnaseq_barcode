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
fastq_file = glob.glob(f"{datadir}*negctrl*.fastq.gz")[0]

#%%
# Read sequences
df_neg = pd.read_csv(
    f"{homedir}/data/extra/neg_ctrl_sequences.tsv", delimiter="\t"
)
# Assign column with sequence length
df_neg = df_neg.assign(seq_len=df_neg["variant"].apply(len))

#%%
# Use skbio to have a generator to iterate over fastq
seqs = skbio.io.read(
    fastq_file,
    format="fastq",
    verify="false",
    variant="illumina1.8",
)

# Set counter
counter = 0

# Define number of samples
n_samples = 10000

print("reading sequences into memory...")
# Initialize list to save sequence objects
seq_list = list()
# Iterate over sequences
for seq in seqs: #itertools.islice(seqs, n_samples):
    if counter % 10000 == 0:
        print(f"reading seq #{counter}")
    # Extract sequence information
    seq_id = seq.metadata["id"]
    sequence = str(skbio.DNA(sequence=seq, validate=False))
    # Append to list
    seq_list.append([seq_id, sequence])
    # Update counter
    counter += 1

# Initialize dataframe to save sequences
names = ["id", "sequence"]
df_seq = pd.DataFrame.from_records(seq_list, columns=names)

# Add sequence length to dataframe
df_seq["seq_len"] = df_seq.sequence.apply(len)

#%%
print("filtering sequences...")
# Define primer sequences
s100_rev = "TACTTTTGATTGCTGTGCCC"
s100_fwd = "ATAACACGGCACGAATAAGC"

# Search for these sequences

# Initialize array to save clone position
fwd_pos = np.empty(len(df_seq))
rev_pos = np.empty(len(df_seq))

# Loop through sequences
for i, seq in df_seq.iterrows():
    # Search clone sequence
    fwd_re = re.search(s100_fwd, seq["sequence"])
    rev_re = re.search(s100_rev, seq["sequence"])
    
    # Save position
    if bool(fwd_re):
        fwd_pos[i] = fwd_re.span()[0]
    else:
        fwd_pos[i] = np.nan
        
    if bool(rev_re):
        rev_pos[i] = rev_re.span()[0]
    else:
        rev_pos[i] = np.nan

# Add column to dataframe
df_seq = df_seq.assign(fwd_prim=fwd_pos, rev_prim=rev_pos)

# Reset index
df_seq.reset_index(inplace=True, drop=True)

#%%
# Filter out sequences without the primer sequences
df_seq = df_seq.dropna()

# Compute distance between primers
df_seq = df_seq.assign(
    primer_dist = df_seq["fwd_prim"] - df_seq["rev_prim"]
)
# Reset index
df_seq.reset_index(inplace=True, drop=True)
#%%
# Filter out sesquences with primer_dist != 170
df_seq = df_seq[df_seq["primer_dist"] == 170]

# Reset index
df_seq.reset_index(inplace=True, drop=True)
#%%

# Define clone binding site
clone = str(
    skbio.DNA("gctagcCAATGCGGgagctc".upper()).reverse_complement()
)

# Initialize array to save clone position
clone_pos = np.zeros(len(df_seq))
# Loop through sequences
for i, seq in df_seq.iterrows():
    # Search clone sequence
    clone_re = re.search(str(clone), seq["sequence"])
    # Save position
    if bool(clone_re):
        clone_pos[i] = clone_re.span()[0]
    else:
        clone_pos[i] = np.nan

# Add column to dataframe
df_seq = df_seq.assign(clone=clone_pos)

# Compute clone distance
clone_dist = (
    (df_seq["clone"] - df_seq["rev_prim"])
    + len(clone)
)

# Select sequences
df_seq = df_seq[(clone_dist == 0) & (clone_pos == 20)]

# Reset index
df_seq.reset_index(inplace=True, drop=True)

#%%

print("Mapping sequence and barcode...")
# Initialize dataframe to save sequences and barcodes
df_neg_map = pd.DataFrame([])

# Extract negative control sequence
df_neg_map["sequence"] = df_seq["sequence"].apply(lambda x: x[60: 60+150])

# Extract barcodes 
df_neg_map["barcode"] = df_seq["sequence"].apply(lambda x: x[0:20])

#%%

print("Mapping sequences to TWIST order...")
# Set boolean array to determine if sequence is listed or not
neg_idx = np.array([False] * len(df_neg_map))

# Loop through each of the sequences
for i, seq in df_neg_map.iterrows():
    # Compute reverse complement of sequence
    s = str(skbio.DNA(seq["sequence"]).reverse_complement())
    # Find if sequence is included in list
    if s in df_neg["variant"].values:
        neg_idx[i] = True

# Filter out sequences that are not exactly as in Guillaume's list
df_neg_map = df_neg_map[neg_idx]

# Reset index
df_neg_map.reset_index(inplace=True, drop=True)

print(f"Number of sequence: {len(df_neg_map)}")

#%%

print("Counting reads for each sequence/barcode pair")

# Generate df with sequence/barcode pairs and their counts
df_counts = df_neg_map.value_counts().reset_index(name="counts")

# Filter out pairs that have <= 3 counts
df_counts = df_counts[df_counts["counts"] >= 3]

#%%
# Write file to memory
print("writing barcode list into memory")
df_counts.to_csv(f"{outputdir}negctrl_barcodes_counts.csv", index=False)


print("Done! Barcodes filtered and quantified")
# %%
