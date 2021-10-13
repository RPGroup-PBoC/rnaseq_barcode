#%%
import os
import glob
import pandas as pd
import git

#%%
# Date for sequencing run of the library to map
DATE = 20210715
# Description to be attached to folder names
DESCRIPTION = '_lacI_titration/'

# Find project parental folder
repo = git.Repo('./', search_parent_directories=True)
homedir = repo.working_dir

# Path to input folder
INPUT_DIR = f'{homedir}/data/demux_sequencing/{DATE}{DESCRIPTION}'

# Path to output folder
OUTPUT_DIR = f'{homedir}/data/processed_sequencing/{DATE}{DESCRIPTION}'

# Generate output directory if it doesn't exist
if not os.path.exists(OUTPUT_DIR):
    os.makedirs(OUTPUT_DIR)

#%%
# Read list of demultiplexed sequences
seq_table = pd.read_csv(INPUT_DIR + 'MANIFEST', comment="#")

# Group sequences by index to feed them into the fastp
seq_group = seq_table.groupby('sample-id')

for group, seq in seq_group:
    print(group)
    # Extract file names
    forward = seq[seq.direction == 'forward'].filename.values[0]
    # reverse = seq[seq.direction == 'reverse'].filename.values[0]
    # Define inputs for fastp
    in1 = f'{INPUT_DIR}{forward}'  # forward input read
    # in2 = f'{INPUT_DIR}{reverse}'  # reverse input read

    # Define outputs
    out1 = f'{OUTPUT_DIR}{forward}'  # reverse input read
    # out2 = f'{OUTPUT_DIR}{reverse}'  # reverse input read
    html_report = f'{OUTPUT_DIR}{DATE}_{group}_fastp_report.html'
    json_report = f'{OUTPUT_DIR}{DATE}_{group}_fastp_report.json'
    report_title = f'{DATE}{DESCRIPTION} fastp report'

    # Define string to be ran on the terminal
    fastp = f'''fastp \
        --in1 {in1} \
        --out1 {out1} \
        --verbose \
        --max_len1 60 \
        --disable_adapter_trimming \
        --html {html_report} \
        --json {json_report} \
        --report_title "{html_report}" \
    '''

    # Run program
    os.system(fastp)
