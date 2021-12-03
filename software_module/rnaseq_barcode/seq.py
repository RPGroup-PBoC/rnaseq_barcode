import pandas as pd
import scipy as sp
import numpy as np
import skbio
import itertools
import regex


def seq2mat(seq,seq_dict):
    mat = sp.zeros((len(seq_dict),len(seq)),dtype=int)
    for i,bp in enumerate(seq):
        mat[seq_dict[bp],i] = 1
    return mat


def choose_dict(dicttype,modeltype='MAT'):
    '''Creates a necessary tool to convert from bp to an index'''
    if dicttype == 'dna':
        seq_dict = {'A':0,'C':1,'G':2,'T':3}
        inv_dict = {0:'A',1:'C',2:'G',3:'T'}
    elif dicttype == 'rna':
        seq_dict = {'A':0,'C':1,'G':2,'U':3}
        inv_dict = {0:'A',1:'C',2:'G',3:'U'}
    elif dicttype == 'protein':
        seq_dict = {
            '*':0,'A':1,'C':2,'D':3,'E':4,'F':5,'G':6,'H':7,'I':8,'K':9,'L':10,
            'M':11,'N':12,'P':13,'Q':14,'R':15,'S':16,'T':17,'V':18,'W':19,'Y':20}
        inv_dict = {v:k for k,v in seq_dict.items()}
    else:
        raise SortSeqError('Unkonwn dicttype: %s'%dicttype)

    if modeltype == 'NBR' or modeltype == 'PAIR':
        seq_dict = {
            ''.join([inv_dict[i],inv_dict[z]]):i*len(seq_dict)+z
            for i in range(len(seq_dict)) for z in range(len(seq_dict))}
        inv_dict = {seq_dict[i]:i for i in seq_dict.keys()}
    return seq_dict,inv_dict


def read_sequences(file, start=0,stop=None):
    """Import sequences from fastq file and return them as dataframe.
    
    Parameters
    ----------
    file : string
        Path to the file containg sequences. Has to be fastq file.
    start : int, default 0
        Index of first sequence to be read.
    stop : int, default None
        Index of last sequence to be read.

    Returns
    -------
    df_seq : DataFrame
        Pandas dataframe containing the sequence as well as sequence length and
        read ID.
    
    """
    # Use skbio to have a generator to iterate over fastq
    seqs = skbio.io.read(
        file,
        format="fastq",
        verify="false",
        variant="illumina1.8"
    )


    # Initialize list to save sequence objects
    seq_list = list()
    
    for seq in itertools.islice(seqs, start, stop):
        # Extract sequence information
        seq_id = seq.metadata["id"]
        sequence = str(skbio.DNA(sequence=seq, validate=False))
        # Append to list
        seq_list.append([seq_id, sequence])

    # Initialize dataframe to save sequences
    names = ["id", "sequence"]
    df_seq = pd.DataFrame.from_records(seq_list, columns=names)

    # Add index and sequence length to dataframe
    df_seq["index"] = file.split("/")[-1]
    df_seq["seq_len"] = df_seq.sequence.apply(len)
    return df_seq


def find_bc(
    seq, 
    ref_seq, 
    ind_filter=None, 
    mismatches=1, 
    barcode_length=5,
    direction="down"
    ):
    """Read barcode from sequence using adjacent reference sequence.
    
    Parameters
    ----------
    seq : string
        Sequence as string.
    ref_seq : string
        Reference sequence used to locate barcode.
    ind_filter : list, default None
        If given, only returns sequences which locate the refernce sequence
        at one of the given positions.
    mismatches : int, default 1
        Number of mismatches allowed in the reference sequence. 
    barcode_length : int, default 5
        Length of the barcode
    direction : string, default "up"
        If "up", the reference sequence is upstream of the barcode in the read.
        If "down", the reference sequence is downstream of the barcode in the read.

    Returns
    -------
    idx : int
        Location of the barcode in the read. If no reference sequence found, or at 
        the wrong location, returns `np.nan`.
    bc : string
        Located barcode. If no reference sequence found, or at 
        the wrong location, returns `np.nan`.
    
    """

    if direction not in ["up", "down"]:
        raise ValueError("`direction` must be either 'up' or 'down'!")

    # Find reference sequence in read
    if mismatches == 0:
        reg = regex.finditer(f"({ref_seq})" + "{s<=0}", seq)
    elif mismatches == 1:
        reg = regex.finditer(f"({ref_seq})" + "{s<=1}", seq)
    else:
        raise ValueError("Only up to one mismatch!")
    # Extract span where reference sequence is located
    match = [m.span() for m in reg]

    # 
    if ind_filter != None and type(ind_filter) not in [list, np.ndarray]:
        ind_filter = [ind_filter]
    
    # Rules for lack of match
    if len(match) == 0:
        idx = np.nan
        bc = np.nan
    else:
        if direction == "down":
            if ind_filter != None:
                if match[0][0] - barcode_length not in ind_filter:
                    idx = np.nan
                    bc = np.nan
                else:
                    idx = match[0][0] - barcode_length
                    bc = seq[match[0][0] - barcode_length : match[0][0]]
            else:
                idx = match[0][0] - barcode_length
                bc = seq[match[0][0] - barcode_length : match[0][0]]
        else:
            if ind_filter != None:
                if match[0][0]  not in ind_filter:
                    idx = np.nan
                    bc = np.nan
                else:
                    idx = match[0][1] + 1
                    bc = seq[match[0][1] : match[0][1] + barcode_length]
            else:
                idx = match[0][1] + 1
                bc = seq[match[0][1] : match[0][1] + barcode_length]

    return idx, bc


def find_close_bc(seq, barcode_list):
    bc_lists = [list(x) for x in barcode_list]
    for j, bc in enumerate(bc_lists):
        if np.sum([x != y for x, y in zip(list(seq), bc)]) == 1:
            return barcode_list[j]
        
    return np.nan

