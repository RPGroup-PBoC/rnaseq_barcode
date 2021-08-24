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
distribution`](https://www.anaconda.com/distribution/) on your computer. 
Instructions on installing Qiime can be found [here](https://docs.qiime2.org/2021.4/install/native/).
Here we show the installation instructions on  for the version of `qiime2` that is used for analysis of data
in this project. There might be newer versions of `qiime2` available, but we cannot guarantee that
the analysis pipeline works for these versions. Note that the only difference between the official instructions
and the ones shown here is the name of the conda environment that is created.
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
This will create a new conda environment and install the necessary libraries to run `qiime2` basic functions.

In general, to run `qiime2` the environment needs to be activated. This is done by writing
in the command line (replace `conda` with `source` in Windows)
```
conda activate qiime2
```

The to switch back the base conda environment, deactivate the current environment by
```
conda deactivate
```

To run our customized python scripts the environment does not need to be activated.
It is done automatically by the `subprocess` module.

The demultiplexing scripts assume the file structure for this project that
looks like:
```
+---code
|   +---processing
|       +---DATE_sequencing_processing_scripts
|           +---demultiplex_seq.py
|
+---data
    +---raw_sequencing
    |    +---DATE_SEQUENCING-INFORMATION
    |        +---sequencing_barcodes_qiime2.tsv (barcodes list)
    |        +---MISEQ (files as given by Illumina)
    |            +---Data
    |                +---Intensities
    |                    +---BaseCalls
    |                        +---R1.fastq.gz (forward read)
    |                        +---R2.fastq.gz (reverse read)
    |                        +---I1.fastq.gz (index read)
```
Where all names written in capital letters are given as variables in the `demultiplex_seq.py` script. The folder structure for
data is quite strict, but can be modified, as long as the `demultiplex_seq.py` is adapted appropriately.

To run the script, simply execute the `demultiplex_seq.py` file from the terminal.

The index barcodes needed for demultplexing are kept in `sequencing_barcodes_qiime2.tsv`. For this file, make sure the
columns are separated by tabs, and not spaces, otherwise `qiime` complains.

## Demultiplexed sequencing processing

After the raw sequences have been split into individual `fastq` files for each
of the indexes, the sequences need to be processed based on their quality,
length, and in the case of paired-end reads, the forward and reverse read must
be stuck together. To perform these tasks, we use
[`fastp`](https://github.com/OpenGene/fastp). The package is installed via
```
conda install -c bioconda fastp
```
After that, all scripts named `processing_seq.py` can be run to process the
short-reads. These scripts assume that the data has been demultiplexed already
as it takes the resulting `fastq` files from the `demux_sequencing` folder
indicated above. The output of these scripts is saved in
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