#%%
import os
import glob
import itertools
import re
import numpy as np
import pandas as pd
import collections
import skbio

#%%
# Define data directory
datadir = '../../../data/processed_sequencing/' +\
          '20190821_operator_library_mapping/'

# Define output dir
outputdir = '../../../data/barcodes/' +\
            '20190821_operator_library_mapping/'

# List all fastq.gz files
fastq_files = glob.glob(f'{datadir}*.fastq.gz')

# Define operator sequences
O1 = skbio.DNA('aattgtgagcggataacaatt'.upper()).reverse_complement()
O2 = skbio.DNA('aaatgtgagcgagtaacaacc'.upper()).reverse_complement()
O3 = skbio.DNA('ggcagtgagcgcaacgcaatt'.upper()).reverse_complement()
operators = {'O1': str(O1),
             'O2': str(O2),
             'O3': str(O3)}

# Define function to map operator from sequence
def op_match(seq):
    '''
    Function to match the operator sequences
    '''
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
        return ['None', 0, 0]

# Define RNAP binding site
rnap = str(skbio.DNA('TTTACACTTTATGCTTCCGGCTCGTATAATGTGTGG').\
           reverse_complement())
# Define sequence post operator
cloning = str(skbio.DNA('GCTAGCGGTGTTTAGTTAGCATCCGGTCTCACATGC').\
          reverse_complement())

def const_match(seq):
    '''
    Function to match the operator sequences
    '''
    # Find RNAP
    rnap_pos = re.search(rnap, seq)
    # Find cloning site
    cloning_pos = re.search(cloning, seq)
    
    # Determine if RNAP was found
    if bool(rnap_pos):
        RNAP = [True] + [*rnap_pos.span()]
    else:
        RNAP = [False, 0, 0]
    # Determine if cloning site was found
    if bool(cloning_pos):
        CLONING = [True] + [*cloning_pos.span()]
    else:
        CLONING = [False, 0, 0]
    
    return RNAP + CLONING
#%%
for i, fastq in enumerate(fastq_files):
    print(fastq)

    # Use skbio to have a generator to iterate over fastq
    seqs = skbio.io.read(datadir + fastq,
                        format='fastq',
                        verify='false',
                        variant='illumina1.8')

    
    print('reading file into memory')
    # Initialize list to save sequences
    seq_list = list()

    # Initialize counter
    counter = 0
    # Iterate over sequences
    for seq in seqs:
        if counter%100000 == 0:
            print(f'read # {counter}')
        # Extract ID
        seq_id = seq.metadata['id']
        # Extract sequence
        seq_str = str(skbio.DNA(sequence=seq, validate=False))
        # Append to list
        seq_list.append([seq_id, f'index{i}', seq_str])      
        counter += 1
    
    # Generate Pandas dataframe from list
    names = ['id', 'index', 'sequence']
    df_seq = pd.DataFrame.from_records(seq_list, 
                                       columns=names)

    # Add length to dataframe
    df_seq['seq_len'] =  df_seq.sequence.apply(len)

    # Map operators
    op_map = list()
    print('mapping operators')
    # Loop through rows
    for seq in df_seq.sequence:
        op_map.append(op_match(seq))

    df_seq = pd.concat([df_seq,
                        pd.DataFrame\
                        .from_records(op_map,
                                      columns=['operator',
                                      'op_begin',
                                      'op_end'])],
                        axis=1)
    print('Done!')

    print('mapping RNAP binding site and cloning site')
    # Initialize list to save results
    const_map = list()
    # Loop through rows
    for seq in df_seq.sequence:
        # Map RNAP and cloning site
        const_map.append(const_match(seq))

    # Generate dataframe from 
    df_const = pd.DataFrame.from_records(const_map,
                                columns=['rnap',
                                'rnap_begin',
                                'rnap_end',
                                'cloning',
                                'cloning_begin',
                                'cloning_end'])

    # Concatenate dataframes
    df_seq = pd.concat([df_seq, df_const], axis=1)
    print('Done!')

    print('writing file into memory')
    # Write file to memory
    df_seq.to_csv(f'{outputdir}index{i+1}_raw_operator_mapping.csv',
                  index=False)

    # Apply strong filter to operators based on length, operator existence and
    # operator position
    print('filtering sequences by length and barcode position')
    df = df_seq[(df_seq.operator != 'None') &
                (df_seq.seq_len == 113) &
                (df_seq.op_begin == 56)]

    # reset index
    df = df.reset_index(drop=True)

    print('filtering sequences by RNAP and cloning site position')
    # Apply filters to sequences based on the presence of the
    # RNAP binding site and the cloning site with the right
    # Position
    df = df[(df.rnap) &
            (df.cloning) &
            (df.rnap_begin == 77) &
            (df.cloning_begin == 20)]

    # Reset index
    df = df.reset_index(drop=True)

    # Group by operator
    df_group = df.groupby('operator')

    # Initialize dataframe to save outcome
    names = ['operator', 'sequence', 'barcode', 'counts']
    df_counts = pd.DataFrame(columns=names)

    # Loop thorough operators
    print('counting unique sequences')
    for group, data in df_group:
        # Count unique barcodes and turn it into a DataFrame
        df_op = data['sequence'].value_counts()\
                                .rename_axis('sequence')\
                                .reset_index(name='counts')
        # Add a column that contains operator
        df_op['operator'] = [group] * len(df_op)
        # Extract barcodes
        df_op['barcode'] = df_op['sequence'].apply(lambda x: x[0:20])
        # Append to dataframe
        df_counts = df_counts.append(df_op,
                                    ignore_index=True,
                                    sort=False)
    print('Done!')

    # Write file to memory
    print('writing barcode list into memory')
    df_counts.to_csv(f'{outputdir}index{i+1}_operator_counts.csv',
                  index=False)
