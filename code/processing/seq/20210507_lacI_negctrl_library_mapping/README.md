---
status: Accepted
reason: Data looks good. 
---

# 2021-05-07 Sequencing run

## Purpose

The purpose of this sequencing run was to map promoters to unique barcodes for three different related (but distinct libraries). The libraries are

1. **Negative Controls**. These are sequences identified by Guillaume as transcriptionally silent. They were ordered in Niko's twist order to be used as negative controls to background subtract residual transcription or DNA contamination from other libraries.
2. **LacUV5 + LacO1 mutants** These promoters have the LacUV5 RNAP site immediately followed by LacO1. This is a pool of several thousand single, double, and higher order mutants in the LacO1 site designed to uniformly sample the energy matrix. This was part of the twist 2020 order and is detailed in [this notebook](https://github.com/RPGroup-PBoC/Reg-Seq2/blob/master/code/experimental_design/twist_order/lacI_titration/generate_sequences.ipynb)
3. **LacUV5 + WT LacO1/2/3** These three promoters are just LacUV5 RNAP with the three Lac operators O1, O2, and O3. These sequences are the same as Manuel and Niko's previous attempt. These promoters were simply ordered as IDT oligos and amplified. 


## Platform
MiSeq system (Thomson lab)

## Sequencing kit
MiSeq **Micro** kit V2. (300 cycles)

## Sequencing Modality
Paired-ends reads, 150 cycles fwd and reverse. Index was also read.

## Materials

| **id** | **barcode-sequence** | **index-primer** | **description** |
| :--: | :--: | :--: | :--: |
| Negative_controls | CACAGT | SC81/GU80 | Transcriptionally silent controls |
| LacO1_mutants | ATGGCT | SC82/GU81 | LacUV5 + LacO1 mutants |
| WT_LacO1/2/3 | CGAGAT | SC83/GU82 | LacUV5 + WT LacO1/2/3 |
