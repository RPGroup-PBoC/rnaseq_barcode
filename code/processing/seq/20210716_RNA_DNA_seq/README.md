---
status: Rejected
reason: Data not clear
---

# 2021-07-15 Sequencing run

## Purpose
In this experiment we sequence the barcodes of a pooled library, which contains
strains of 7 different repressor copy numbers (LacI), each with in individual barcode, with
and additional negative control strain. Each strain is a library of three operators for
the repressor LacI. A random barcode maps to the operator. In a previous experiment, we
mapped the random barcode to the operator, such that we only need to sequence the
barcode to identify the operator.

All files used for analysis of this dataset are included in this folder.

## Platform
MiSeq system (Thomson lab)

## Sequencing kit
MiSeq Kit V2.  (50 cycles)

## Sequencing Modality
Single-end read (50 cycles Read 1, 8 cycles Index Read 1)

## Materials

| **id** | **barcode-sequence** | **biological replicate** | **technical replicate** | **type** |
| :--: | :--: | :--: | :--: | :--: |
| 1_1_DNA | CACAGT | 1 | 1 | gDNA |
| 1_2_DNA |	ATGGCT | 1 | 1 | gDNA | 
| 2_1_DNA |	CGAGAT | 2 | 2 | gDNA |
| 2_2_DNA |	ACACTG | 2 | 2 | gDNA |
| 1_1_RNA |	CATTCG | 1 | 1 | cDNA |
| 1_2_RNA |	GCATAG | 1 | 1 | cDNA |
| 2_1_RNA |	ACTAGC | 2 | 2 | cDNA |
| 2_2_RNA |	CAGTAC | 2 | 2 | cDNA |


| **strain** | **gfp-barcode-sequence** | **reverse_complement_gfp_barcode_sequence** | **repressor** |
| :--: | :--: | :--: | :--: | :--: |
| Negative Controls | GTGTT | AACAC | 0 |
| R0 | ACCTT | AAGGT | 0 |
| R22 | TCAGT | ACTGA | 22 |
| R60 | GGAAT | ATTCC | 60 |
| R124 | ATTGG | CCAAT | 124 |
| R260 | AAGCG | CGCTT | 260 |
| R1220 | TTCTC | GAGAA | 1220 |
| R1740 | ATGAC | GTCAT | 1740 |