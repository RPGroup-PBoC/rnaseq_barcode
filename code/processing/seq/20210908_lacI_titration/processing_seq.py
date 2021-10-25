#%%
import os
import glob
import pandas as pd
import git
import multiprocessing
#%%
# Date for sequencing run of the library to map
DATE = 20210908
# Description to be attached to folder names
DESCRIPTION = '_lacI_titration/'

# Find project parental folder
repo = git.Repo('./', search_parent_directories=True)
homedir = repo.working_dir

# Path to input folder
datadir = homedir + f'/data/raw_sequencing/{DATE}{DESCRIPTION}'

# Find sequencing files
files = glob.glob(datadir + '/Project_*/Sample_*/*.fastq.gz')
# Generate output directory if it doesn't exist
OUTPUT_DIR = homedir + f'/data/processed_sequencing/{DATE}{DESCRIPTION}'
if not os.path.exists(OUTPUT_DIR):
    os.makedirs(OUTPUT_DIR)

# Group sequences by index to feed them into the fastp


def run_fastp(file):
    # Define inputs for fastp
    in1 = file  # forward input read
    name = file.split("/")[-1]
    print(in1)
    # Define outputs
    out1 = f'{OUTPUT_DIR}{name}'
    print(out1)
    html_report = f'{OUTPUT_DIR}{DATE}_{name}_fastp_report.html'
    json_report = f'{OUTPUT_DIR}{DATE}_{name}_fastp_report.json'
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


a_pool = multiprocessing.Pool(40)
result = a_pool.map(run_fastp, files)