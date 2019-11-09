import os
import numpy as np
import pandas as pd
import glob
import imp
import git
import rnaseq_barcode.flow as flow

# Set the experiment constants from folder name
dirname = os.getcwd().split('/')[-1]
DATE = int(dirname.split('_')[0])
RUN_NO = int(dirname.split('_')[1][1:])
USERNAME = 'mrazomej'
gating_fraction = 0.4

# Find project parental directory
repo = git.Repo('./', search_parent_directories=True)
homedir = repo.working_dir
# Load all files.
files = glob.glob(f'{homedir}/data/flow/csv/{DATE}*_r{RUN_NO}*.csv')

# Set up the DataFrame
colnames = ['date', 'username', 'phenotype', 'operator', 'strain', 'IPTGuM',
            'mean_FITC_H']
df = pd.DataFrame([], columns=colnames)

for f in files:
    print(f)
    # Extra feature for lost data, replaced by the 0.1 µM IPTG (See README.md)
    if '20191105_r1_wt_O3_auto_0uMIPTG.csv' in f:
        fnew = f.replace('0uM', '0.1uM')
        print(f'data replaced for {fnew}')
        # Load in the data
        data = pd.read_csv(fnew)

    elif '20191105_r1_wt_O3_delta_0uMIPTG.csv' in f:
        fnew = f.replace('0uM', '0.1uM')
        print(f'data replaced for {fnew}')
        # Load in the data
        data = pd.read_csv(fnew)

    else:
        # Load in the data
        data = pd.read_csv(f)

    # Get the identifying finformation.
    date, _, phenotype, operator, strain, conc = f.split('/')[-1].split('_')

    conc = float(conc.split('uM')[0])
    if (strain != 'auto') and (strain != 'delta'):
        rep = int(strain.split('R')[-1])
    else:
        rep = 0

    # Gate data
    gated = flow.gaussian_gate(data, gating_fraction)

    # Compute the mean
    mean_FITC = gated['FITC-H'].mean()

    # Assemble the dictionary
    samp_dict = dict(date=date, username=USERNAME, phenotype=phenotype,
                     operator=operator, strain=strain, IPTGuM=conc,
                     mean_FITC_H=mean_FITC, repressors=rep)
    df = df.append(samp_dict, ignore_index=True)


fc_dfs = []
grouped = df[df['strain']!='auto'].groupby(['IPTGuM', 'operator'])
mean_auto_df = df[df['strain'] == 'auto']
for g, d in grouped:
    mean_auto = mean_auto_df[mean_auto_df['IPTGuM']==g[0]]['mean_FITC_H'].values[0]
    mean_delta = d.loc[d['strain'] == 'delta']['mean_FITC_H'].values[0]
    d['fold_change'] = (d['mean_FITC_H'] - mean_auto) /  (mean_delta - mean_auto)
    fc_dfs.append(d)

fold_change_df = pd.concat(fc_dfs, axis=0)

if not os.path.exists('./output/'):
    os.mkdir('./output/')
# Save to a CSV.
fold_change_df.to_csv(f'output/{DATE}_r{RUN_NO}_fold_change.csv', index=False)
