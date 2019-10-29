## `templates`
In this directory, we store templates of experiment descriptions and processing
scripts such that updating them for each experiment is simple.

## Raw sequencing processing

- `demultiplex_seq.py` : This script takes raw `fastq` files from the
  sequencing machine and splits them into different files according to a list
  of sequencing index. See `processing/README.md` for more information.
- `sequencing_barcodes_qiime2.tsv` : Example tab separated value file listing
  the sequencing index nucleotides to split files. This specific format is
  required by the `qiime2` package.
- `processing_seq.py` : This script runs the `fastp` processing pipeline to
  filter the sequencing data, merge paired-end reads, among other things. See
  `processing/README.md` for more information.
- `sequencing_README.md` : markdown file to input information about a
  sequencing experiment.


## Raw flow cytometry processing

- `flow_README.md` : markdown file o input information about flow cytometry
  experiment.
- `YYYYMMDD_rx_fcsrename.csv` : Standard `csv` file used by the shell script
  `fcs_rename.py` to rename `fcs` files in bulk. Please make sure to use the
  standardize filename format:
  ```
  YYYYMMDD_run#_phenotype_operator_rbs_XuMIPTG
  ```
  This file must have the same spatial arrangement that the 96 well plate that
  went into the cytometer had for the files to be properly renamed.
- `processing_flow.py` : Script that takes the `csv` parsed flow cytometry data
  and computes the mean expression level along with the fold-change in gene
  expression given an autofluorescence strain and a ∆*lacI* strain.
- `analysis_flow.py` : Script that uses the output from `processing_flow.py`
  and displays the experimental fold-change along with the theoretical
  predictions to assess the quality of the data.