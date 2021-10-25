---
status: Rejected
reason: Data too noisy
---

# 2021-09-08 Sequencing run

## Purpose
This experiment was a follow-up experiment to the sequencing run from July 15th.
In that sequencing run we got noisy data and expected that we did not have
enough reads to get the desired resolution. Hence, we used the same library,
and had to repeat the second PCR amplification. This time we only used the second
technical replicate of each biological replicate, and used the sequencing center at
Caltech, which has a HiSeq available, which can give up to 300 million reads.

## Platform
HiSeq system (Caltech)

## Sequencing Modality
Single-end read (75 cycles Read 1, 8 cycles Index Read 1)

## Materials

| **id** | **barcode-sequence** | **biological replicate** | **technical replicate** | **type** |
| :--: | :--: | :--: | :--: | :--: |

| 1_2_DNA |	ATGGCT | 1 | 1 | gDNA | 
| 2_2_DNA |	ACACTG | 2 | 2 | gDNA |
| 1_2_RNA |	GCATAG | 1 | 1 | cDNA |
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

