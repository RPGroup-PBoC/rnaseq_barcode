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
fastq_file = glob.glob(f"{datadir}*WTlac*.fastq.gz")[0]

#%%
# Define operator sequences
# Forward operators
O1_rev = skbio.DNA("aattgtgagcggataacaatt".upper())
O2_rev = skbio.DNA("aaatgtgagcgagtaacaacc".upper())
O3_rev = skbio.DNA("ggcagtgagcgcaacgcaatt".upper())
# Reverse complement
O1 = O1_rev.reverse_complement()
O2 = O2_rev.reverse_complement()
O3 = O3_rev.reverse_complement()
operators = {
    "O1": str(O1),
    "O2": str(O2),
    "O3": str(O3),
    "O1_rev": str(O1_rev),
    "O2_rev": str(O2_rev),
    "O3_rev": str(O3_rev),
}

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


def op_match(seq):
    """
    Function to match the operator sequences
    """
    # Loop through operators
    for key, item in operators.items():
        # Find operator and return boolean if found
        op_pos = re.search(item, seq)
        # If found return the operator and break loop
        if bool(op_pos):
            return [key] + [*op_pos.span()]
            break

    # If none match, return none
    if not bool(op_pos):
        return ["None", 0, 0]


def rev_comp(seq):
    """
    Function that takes a string, converts it into skbio.DNA
    and takes the reverse complement
    """
    return str(skbio.DNA(seq, validate=False).reverse_complement())


#%%
print("Reading WTlac sequences into memory")

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
for seq in seqs:  # itertools.islice(seqs, n_samples):
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

print("Mapping operator sequences")
# Map operators
op_map = list()
# Loop through rows
for seq in df_seq.sequence:
    op_map.append(op_match(seq))

df_seq = pd.concat(
    [
        df_seq,
        pd.DataFrame.from_records(
            op_map, columns=["operator", "op_begin", "op_end"]
        ),
    ],
    axis=1,
)

# Reverse complement sequences which had a reversed operator
# Find forward sequences
bool_forward = ["_rev" in x and x != "None" for x in df_seq.operator]

# Reverse complement forward sequences
df_seq.loc[bool_forward, "sequence"] = [
    rev_comp(seq) for seq in df_seq[bool_forward]["sequence"]
]

# Remap operators after having reversed sequences
op_map = list()
# Loop through rows
for seq in df_seq.sequence:
    op_map.append(op_match(seq))

df_seq[["operator", "op_begin", "op_end"]] = pd.DataFrame.from_records(
    op_map, columns=["operator", "op_begin", "op_end"]
)
print("Done mapping operators...")

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
print("2. Filtering by operator sequence and position")

# Copy filtered step
df = df_filt

# Remove sequences with no mapped operator
df_filt = df[df["operator"] != "None"]
# Store rejected sequences
df_r = df[df["operator"] == "None"][["id", "sequence"]]
df_r = df_r.assign(reject_by="operator_seq")
df_reject = df_reject.append(df_r, ignore_index=True)

# Reset index
df_filt.reset_index(inplace=True, drop=True)
df_reject.reset_index(inplace=True, drop=True)

# Compute the distance between the end of the primer and the operator
op_dist = (df_filt["op_begin"] - df_filt["rev_prim"]) - len(rev_prim)

# Copy filtered step
df = df_filt

# Select sequences
df_filt = df[op_dist == 0]
# Store rejected sequences
df_r = df[op_dist != 0][["id", "sequence"]]
df_r = df_r.assign(reject_by="operator_pos")
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
print("3. Filtering by RNAP binding site sequence and position")

# Initialize array to save RNAP position
rnap_pos = np.zeros(len(df_filt))
# Loop through sequences
for i, seq in df_filt.iterrows():
    # Search RNAP sequence
    rnap_re = re.search(str(rnap), seq["sequence"])
    # Save position
    if bool(rnap_re):
        rnap_pos[i] = rnap_re.span()[0]

# Add column to dataframe
df_filt = df_filt.assign(rnap=rnap_pos)

# Compute RNAP distance
rnap_dist = (
    (df_filt["rnap"] - df_filt["rev_prim"])
    - len(rev_prim)
    - len(operators["O1"])
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
print("4. Filtering by cloning site sequence and position")

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

# Group by operator
df_group = df_filt.groupby("operator")

# Initialize dataframe to save outcome
names = ["operator", "sequence", "barcode", "counts"]
df_counts = pd.DataFrame(columns=names)

# Loop thorough operators
for group, data in df_group:
    # Count unique barcodes and turn it into a DataFrame
    df_op = (
        data["sequence"]
        .value_counts()
        .rename_axis("sequence")
        .reset_index(name="counts")
    )
    # Add a column that contains operator
    df_op = df_op.assign(operator=group)
    # Extract barcodes
    df_op["barcode"] = df_op["sequence"].apply(lambda x: x[0:20])
    # Append to dataframe
    df_counts = df_counts.append(df_op, ignore_index=True, sort=False)

# Write file to memory
print("writing barcode list into memory")
df_counts.to_csv(f"{outputdir}WTlac_barcodes_counts.csv", index=False)

df_reject.to_csv(f"{outputdir}WTlac_rejected_seq.csv", index=False)

print("Done! Barcodes filtered and quantified")
# %%
