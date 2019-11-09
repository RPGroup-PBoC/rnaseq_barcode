#%%
import os
import glob
import pandas as pd
import git

#%%
# Date for sequencing run of the library to map
DATE = 20190821
# Description to be attached to folder names
DESCRIPTION = '_operator_library_mapping/'

# Path to input folder
# Find project parental directory
repo = git.Repo('./', search_parent_directories=True)
homedir = repo.working_dir
INPUT_DIR = f'{homedir}/data/demux_sequencing/{DATE}{DESCRIPTION}'

# Path to output folder
OUTPUT_DIR = f'{homedir}/data/processed_sequencing/{DATE}{DESCRIPTION}'

# Generate output directory if it doesn't exist
if not os.path.exists(OUTPUT_DIR):
    os.makedirs(OUTPUT_DIR)

#%%
# Read list of demultiplexed sequences
seq_table = pd.read_csv(INPUT_DIR + 'MANIFEST')

# Group sequences by index to feed them into the fastp
seq_group = seq_table.groupby('sample-id')

for group, seq in seq_group:
    # Extract file names
    forward = seq[seq.direction == 'forward'].filename.values[0]
    reverse = seq[seq.direction == 'reverse'].filename.values[0]
    # Define inputs for fastp
    in1 = f'{INPUT_DIR}{forward}'  # forward input read
    in2 = f'{INPUT_DIR}{reverse}'  # reverse input read

    # Define outputs
    merged_out = f'{OUTPUT_DIR}{group}_merged.fastq.gz'  # merged reads
    html_report = f'{OUTPUT_DIR}{DATE}_{group}_fastp_report.html'
    json_report = f'{OUTPUT_DIR}{DATE}_{group}_fastp_report.json'
    report_title = f'{DATE}{DESCRIPTION} fastp report'

    # Define string to be ran on the terminal
    fastp = f'''fastp \
        --in1 {in1} \
        --in2 {in2} \
        --merge \
        --merged_out {merged_out} \
        --include_unmerged \
        --verbose \
        --disable_length_filtering \
        --correction \
        --html {html_report} \
        --json {json_report} \
        --report_title "{html_report}" \
        --max_len1 113 \
        --max_len2 113 \
        --overlap_len_require 103 \
        --overlap_diff_percent_limit 5
    '''

    # Run program
    os.system(fastp)
