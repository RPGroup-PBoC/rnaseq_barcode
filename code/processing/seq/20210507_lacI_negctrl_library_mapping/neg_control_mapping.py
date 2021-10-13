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
fastq_file = glob.glob(f"{datadir}*Negative_controls*.fastq.gz")[0]

#%%

# Forward primer
fwd_prim = skbio.DNA("GCTTATTCGTGCCGTGTTAT").reverse_complement()
# Reverse primer
rev_prim = skbio.DNA("GGGCACAGCAATCAAAAGTA").reverse_complement()

# Define clone binding site
clone = str(skbio.DNA("gctagcCAATGCGGgagctc".upper()).reverse_complement())

#%%
print("Reading negative control sequences into memory")

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
for seq in seqs:#itertools.islice(seqs, n_samples):
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

df_cont = pd.read_csv(f"{homedir}/data/extra/Negative_Control_Promoters_best.txt", delimiter="\t")
neg_controls = [str(skbio.DNA(x.upper())) for x in df_cont.variant]
neg_controls_rev = [str(skbio.DNA(x.upper()).reverse_complement()) for x in df_cont.variant]
controls = {}
controls_rev = {}
for i, x in enumerate(neg_controls_rev):
    controls[i] = x
    controls_rev[x] = i

controls_fwd = {}
for i, x in enumerate(neg_controls):
    controls_fwd[x] = i


def rev_comp(seq):
    """
    Function that takes a string, converts it into skbio.DNA
    and takes the reverse complement
    """
    return str(skbio.DNA(seq, validate=False).reverse_complement())


def op_match(seq):
    """
    Function to match the control sequences
    """
    # Loop through control
    for key, item in controls.items():
        # Find control and return boolean if found
        op_pos = re.search(item, seq)
        # If found return the control and break loop
        if bool(op_pos):
            return [key] + [*op_pos.span()]
            break

    # If none match, return none
    if not bool(op_pos):
        return ["None", 0, 0]
    
print("Mapping control sequences")
# Map control
op_map = list()
# Loop through rows

for seq in df_seq.sequence:
    op_map.append(op_match(seq))

df_seq = pd.concat(
    [
        df_seq,
        pd.DataFrame.from_records(
            op_map, columns=["control", "cont_begin", "cont_end"]
        ),
    ],
    axis=1,
)

# Reverse complement sequences which had a reversed operator
# Find forward sequences
bool_forward = [x != "None" for x in df_seq.control]

# Reverse complement forward sequences
df_seq = df_seq[bool_forward]
#df_seq['sequence'] = [rev_comp(seq) for seq in df_seq["sequence"]]
df_seq.reset_index(inplace=True)

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


# Filtering sequences
print("Filtering sequences by:")

print("1. Filtering by separation between primers")

# Save original dataframe length
len_df = len(df_seq)

# Filter by length
df_filt = df_seq[df_seq["prim_dist"] == 170]

# Save rejected sequences
df_reject = df_seq[df_seq["prim_dist"] != 170][["id", "sequence"]]
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


# Reset index
df_filt.reset_index(inplace=True, drop=True)

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


# Extract barcodes and control
df_op = df_op.assign(
    barcode=df_op["sequence"].apply(lambda x: x[0:20]),
    control=df_op["sequence"].apply(lambda x: x[60:210]),
)
print(df_op)
df_op = df_op.assign(
    control_id=[controls_rev[x] if x in controls_rev.keys() else "none" for x in df_op.control]
)
# Print percentage of sequences removed
print(
    f"""
Cumulative percentage of original sequences removed: 
{np.round(100 - len(df_op) / len_df * 100, 2)}%
"""
)
df_op = df_op[df_op.control_id != 'none']
# Write file to memory
print("writing barcode list into memory")
df_op.to_csv(f"{outputdir}neg_controls_barcodes_counts.csv", index=False)

df_reject.to_csv(f"{outputdir}neg_controls_rejected_seq.csv", index=False)

print("Done! Barcodes filtered and quantified")
# %%
