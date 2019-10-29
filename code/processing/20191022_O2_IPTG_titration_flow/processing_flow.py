import numpy as np
import pandas as pd
import glob
import imp
import rnaseq_barcode.flow as flow

# Define the experiment parameters
DATE = 20191022
RUN_NO = 1
USERNAME = 'mrazomej'
gating_fraction = 0.4

# Load all files.
files = glob.glob(f'../../../data/flow/csv/{DATE}*_r{RUN_NO}*.csv')

# Set up the DataFrame
colnames = ['date', 'username', 'phenotype', 'operator', 'strain', 'IPTGuM',
            'mean_FITC_H']
df = pd.DataFrame([], columns=colnames)

for f in files:
    print(f)
    # Get the identifying finformation.
    date, _, phenotype, operator, strain, conc = f.split('/')[-1].split('_')
    conc = float(conc.split('uM')[0])
    if (strain != 'auto') and (strain != 'delta'):
        rep = int(strain.split('R')[-1])
    else:
        rep = 0
    # Load in the data
    data = pd.read_csv(f)
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

# Save to a CSV.
fold_change_df.to_csv(f'output/{DATE}_r{RUN_NO}_fold_change.csv', index=False)
