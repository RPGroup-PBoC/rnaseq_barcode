# `processing`

This folder contains all code executed to transform or generate data. Within
this directory, some subdirectories contain summaries of specific experiments,
the code used to process and transform the data, and when indicated the output
of any processing functions.

# Sequencing data
The following steps are taken to transform raw `fastq` files into tidy data
formats that can be easily handled in Python for further analysis.

## Sequencing demultiplexing

All scripts named `demultiplex_seq.py` split `fastq` files given a series of
index reads. To demultiplex our sequencing runs we used the
[`qiime2`](https://qiime2.org) platform. The `qiime2` developing team strongly
suggests installing the platform in a separate `conda` environment. To do so,
there must be already a functional [`Anaconda
distribution`](https://www.anaconda.com/distribution/) on your computer. Then to
generate the environment (named `qiime2` in our scripts, but that is easy to
change), you must run

**macOS/OS X (64-bit)**
```
conda update conda
conda install wget

wget https://data.qiime2.org/distro/core/qiime2-2019.7-py36-osx-conda.yml
conda env create -n qiime2 --file qiime2-2019.7-py36-osx-conda.yml

rm qiime2-2019.7-py36-osx-conda.yml
```

**Linux (64-bit)**
```
conda update conda

wget https://data.qiime2.org/distro/core/qiime2-2019.7-py36-linux-conda.yml
conda env create -n qiime2 --file qiime2-2019.7-py36-linux-conda.yml

rm qiime2-2019.7-py36-linux-conda.yml
```
This will install the necessary libraries to run `qiime2` basic functions.

Normally to run `qiime2` you would have to activate the environment by writing
in the command line
```
source activate qiime2
```
But our customized python scripts do this for you internally using the
`subprocess` module. So there is no need to activate the environment before
running our sequencing demultiplexing scripts.

The demultiplexing scripts assume the file structure for this project that
looks like:
```
+---code
|   +---processing
|       +---DATE_sequencing_processing_scripts
|           +---DATE_demultiplex_seq.py
|
+---data
    +---raw_sequencing
    |    +---DATE_sequencing_data
    |        +---sequencing_barcodes_qiime2.tsv (barcodes list)
    |        +---miseq_output (files as given by Illumina)
    |            +---Data
    |                +---Intensities
    |                    +---BaseCalls
    |                        +---R1.fastq.gz (forward read)
    |                        +---R2.fastq.gz (reverse read)
    |                        +---I1.fastq.gz (index read)
    +---demux_sequencing
        +---DATE_demultiplex_data (demultiplex fastq here)
            +---qiime2_output (output objects form qiime2)
            +---tmp
```

## Demultiplexed sequencing processing

After the raw sequences have been split into individual `fastq` files for each
of the indexes, the sequences need to be processed based on their quality,
length, and in the case of paired-end reads, the forward and reverse read must
be stuck together. To perform these tasks, we use the recently published tool
[`fastp`](https://github.com/OpenGene/fastp). The installation of this tool,
just as for `qiime2` can be done via `conda`. All that needs to be done is type
in the terminal
```
conda install -c bioconda fastp
```
After that, all scripts named `processing_seq.py` can be ran to process the
short-reads. These scripts assume that the data has been demultiplexed already
as it takes the resulting `fastq` files from the `demux_sequencing` folder
indicated above. The output of these scripts is saved under
```
+---data
    +---processed_sequencing
        +---DATE_experiment_description
```
Specifically `fastp` generates individual `fastq` files with all the reads
(merged for paired-end runs) for each index. It also generates `HTML` and
`JSON` summaries of the sequencing processing, listing average read quality,
base composition per position, among other useful quantities.


# Flow Cytometry processing

As control experiments for our input-output functions, we measured known
strains' gene expression with a flow cytometer. The raw output from the
cytometer are `fcs` files, a standard format for this technology. To
parse these files and convert them into `csv` files we use the Python library
`fcsparser`. To install such library with **Anaconda** run the following command
in the terminal
```
conda install -c bioconda fcsparser
```
With this library installed then the user can use the `fcs_processing.py`
script in `code/shell_scripts/` to export the data from `fcs` to `csv`.

## Automatic gating and fold-change computation

A 2D Gaussian is fit to the front and side scattering data to gate the flow data
using an unsupervised method. From there, a fraction of 40% of the data is kept
assuming that everything that falls outside of this range must be either not
single cells or other debris.

The script `processing_flow.py` performs this automatic gating and computes the
mean fluorescence signal to then compute the fold-change in gene expression
given an autofluorescence and a ∆*lacI* strain. This script assumes that all of
the `csv` files containing the flow-cytometer data are contained in
```
+---data
    +---flow
        +---csv
```
And are only distinguished by the filenames containing the date and the run
number of the experiment. The file naming format for both the `csv` and the
`fcs` files is standardized as
 ```
 YYYYMMDD_run#_phenotype_operator_rbs_XuMIPTG
 ```
The script can easily parse the information from each measurement and add it
into a tidy dataframe. In addition, this script has as an output a `csv` file
containing computed fold-changes for all strains in the experiment.

After the data is processed, the script `analysis_flow.py` plots the data next
to the theoretical predictions. This serves as an assessment of the quality of
the data.