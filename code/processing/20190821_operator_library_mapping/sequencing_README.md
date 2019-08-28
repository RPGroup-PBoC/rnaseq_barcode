---
status: Rejected
reason: Experiment not yet completed 
---

# 2019-08-21 Sequencing run

## Purpose
Mapping the O1, O2 and O3 oprator library barcodes to the corresponding set 
of random barcodes.

## Platform
MiSeq system (Thomson lab)

## Sequencing kit
MiSeq nano kit V2. (300 cycles)

## Sequencing Modality
Paired-ends reads

## Materials

| **id** | **barcode-sequence** | **description** |
| :--: | :--: | :--: |
| index-1 | CAGTACAT | Technical replica for sequencing |
| index-2 | ACTAGCAT | Technical replica for sequencing |

## Notes and Observations
This is our very first sequencing run. We were helped by Jeff Park from the 
Thomson lab. He taught us how to operate the MiSeq and how to load the samples
into the machine


## Library Preparation

1. The general idea behind this barcode mapping run is that we have different
   operators, each of which is to be "linked" to a barcode (or many barcodes).
   The first step is to order gBLOCKs from IDT that encode our three wildtype
   operator sequences.
2. The three constructs ordered are as follows: 
    - O1 Construct: 
      `ACCTGTAATTCCAAGCGTCTCGAGTTTACACTTTATGCTTCCGGCTCGTATAATGTGTGGAATTGTGAGCGGATAACAATTGCTAGCGGTGTTTAGTTAGCATCCGGTCTCACATGCTAGGCAAC`
    - O2 Construct: 
    `ACCTGTAATTCCAAGCGTCTCGAGTTTACACTTTATGCTTCCGGCTCGTATAATGTGTGGAAATGTGAGCGAGTAACAACCGCTAGCGGTGTTTAGTTAGCATCCGGTCTCACATGCTAGGCAAC`
    - O3 Construct: `ACCTGTAATTCCAAGCGTCTCGAGTTTACACTTTATGCTTCCGGCTCGTATAATGTGTGGGGCAGTGAGCGCAACGCAATTGCTAGCGGTGTTTAGTTAGCATCCGGTCTCACATGCTAGGCAAC`
3. After ordering the gBLOCKs, the next step is to add the barcodes via
   overhang PCR reactions. The number of cycles to be used for this barcoding
   PCR reaction were determined by quantitative PCR with KAPA SYBR Fast; cycles
   were selected at exponential growth in the qPCR data so as to not saturate
   the barcodes and introduce bias.
4. After checking the number of cycles to be used for barcoding via qPCR, we
   next performed the actual PCR with the selected number of cycles. Bands were
   visualized on a 1% agarose gel and PCR reactions purified with a Zymo Clean
   and Concentrate kit.
5. Primers for this PCR reaction also add restriction site overhangs; in the
   next step, pLib vector was digested with SalI and SbfI and the PCR products
   were digested with SbfI and XhoI. After digestion, pLib and each digested
   PCR product were purified again with the Zymo Clean and Concentrate kit.
6. A T7 ligation was performed with the digested pLib and each digested
   Barcode-Operator PCR product, and incubated for 1 hour at room temperature.
   The three ligation mixtures were then purified with the Zymo Clean and
   Concentrate kit.
7. NEB 2989K electrocompetent cells were then used for electroporation; 1ul of
   each cleaned ligation product was transformed. After each transformation,
   cells were plated, then diluted 100x and plated again (6 plates in total,
   oversized agar pads with kanamycin at 25µg/mL).
8. Plates were incubated at 30 degrees Celsius for 24 hours.
9. Colonies from each plate were scraped into individual tubes and glycerol
   stocked. 5 million colonies of the "scrapings" for each of the three
   different operators were also pooled into 400mL LB, grown for 16 hours, and
   then the Barcode-Operator-pLib plasmids maxiprepped and then subsequently
   cleaned with a Zymo Clean and Concentrate kit.
10. After collecting the purified pLib vectors, we proceeded with preparations
    for the NGS DNA sequencing.

### Library preparation for NGS MiSeq
11. qPCR was performed on the pooled, 100x diluted barcoded-operator libraries
    with primers GU60 and GU79 and 1ng of DNA, with Luna qPCR 2x Master Mix. 10
    cycles was selected for subsequent PCR.
12. Standard PCR with Q5 Master Mix was performed on the 100x diluted
    barcoded-operator libraries with primers GU60 and GU79 for 12 cycles with
    1ng DNA.
13. 1µL of PCR product was run on an agarose gel to verify band size. PCR
    reactions were purified with a Zymo Clean and Concentrate kit.
14. A second qPCR was then performed to check the amplification of the
    index-adding primers. This was performed on the PCR amplicons from the
    previous step. Two indexes were tested: GU70/GU86 and GU70/GU87. 8 cycles
    appeared to be optimal.
15. PCR was performed with 8 cycles and the two sets of primers (GU70/GU86 and
    GU70/GU87), using 1ng of PCR amplicons from Step 12 as templates. 
16. 1µL of PCR product was run on an agarose gel to verify band size. PCR
    reactions were purified with a Zymo Clean and Concentrate kit.
17. Purified, indexed libraries were quantified on a Qubit fluorometer and
    normalized according to Illumina MiSeq standard 'DNA Library Preparation'
    instructions.