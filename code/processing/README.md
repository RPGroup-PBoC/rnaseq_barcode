# `processing`

This folder contains all code executed to transform or generate data. Within
this directory, there are subdirectories that contain summaries of specific
experiments, the code used to process and transform the data, and when
indicated the output of any processing functions.

## Sequencing demultiplexing

All scripts named `DATE_demultiplex_seq.py` split `fastq` files given a
series of index reads. To demultiplex our sequencing runs we used the
[`qiime2`](https://qiime2.org) platform. The `qiime2` developing team strongly
suggests installing the platform on a separate `conda` environment. To do so
there must be already a functional [`Anaconda
distribution`](https://www.anaconda.com/distribution/) in your computer. Then
to generate the environment (named `qiime2` in our scripts, but that is easy to
change) you must run

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