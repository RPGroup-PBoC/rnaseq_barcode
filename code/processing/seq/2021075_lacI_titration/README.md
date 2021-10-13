---
status: Rejected
reason: Experiment not yet completed 
---

# 2019-11-08 Sequencing run

## Purpose
This is an experiment, conducted on an Illumina MiSeq in Matt Thomson's
laboratory. The run set-up was assisted by Jeff Park. We used an Illumina MiSeq
v2 kit (300 cycles), which enables up to 15 million reads. The purpose of this
experiment, which is a pooled RNA-seq/DNA-seq experiment, is to determine
whether the "simple repression motif", which the Phillips lab has previously
studied using flow cytometry, fluorescence microscopy, single-molecule RNA
FISH, and other methods, can be reproduced using sequencing. RNAseq reads will
be normalized to DNAseq reads for each sample, and these values will be
"mapped" onto plots made by Hernan Garcia and Rob Phillips (2011 PNAS). A 20N
barcode determines which operator is present (O1, O2, or O3) and at least 1000
barcodes were present for each operator after the initial "mapping" run. An
additional 4-nt "GFP barcode" is used to map the strain in which the operator
libraries were cloned.

## Platform
MiSeq system (Thomson lab)

## Sequencing kit
MiSeq Kit V2. (300 cycles)

## Sequencing Modality
Paired-end read (151 cycles Read 1, 8 cycles Index Read 1, 50 cycles Read 2)

## Materials

| **id** | **barcode-sequence** | **biological replica** | **technical replica** | **type** |
| :--: | :--: | :--: | :--: | :--: |
| B01_R01_cDNA | CACAGTAT | 01 | 01 | cDNA |
| B01_R02_cDNA | ATGGCTAT | 01 | 02 | cDNA |
| B02_R01_cDNA | CGAGATAT | 02 | 01 | cDNA |
| B02_R02_cDNA | ACACTGAT | 02 | 02 | cDNA |
| B01_R01_gDNA | CATTCGAT | 01 | 01 | gDNA |
| B01_R02_gDNA | GCATAGAT | 01 | 02 | gDNA |
| B02_R01_gDNA | ACTAGCAT | 02 | 01 | gDNA |
| B02_R02_gDNA | CAGTACAT | 02 | 02 | gDNA |

| **strain** | **gfp-barcode-sequence** | **reverse_complement_gfp_barcode_sequence** | **barcode_number** | **repressor** |
| :--: | :--: | :--: | :--: | :--: |
| hg105 | ACGT | ACGT | 09 | R0 |
| hg104 | AAGC | GCTT | 10 | R22 |
| rbs1147 | AGTC | GACT | 12 | R60 |
| rbs446 | TCAG | CTGA | 18 | R124 |
| rbs1027 | GTAC | GTAC | 19 | R260 |
| rbs1 | TGCA | TGCA | 20 | R1220 |
| 1l | CGTT | AACG | 21 | R1740 |

## Notes and Observations
Insert details on analysis pipeline here (e.g. using only 60 reads, etc)

## Library Preparation
A RiboJ::sfGFP sequence (obtained from the Sri Kosuri lab, UCLA) was PCR
amplified with NM08 (Fwd primer) and NM09, NM10, NM12, NM18, NM19, NM20 or NM21
(Reverse primers, which add on the 4-nt barcodes unique to each final strain). 

The operator-barcode libraries (sequenced in the last experiment) were digested
with NheI and BsaI. The PCR-amplified riboJ::sfGFP constructs were digested
with BsaI and NcoI, and then ligated into the digested pLib operator-barcode
libraries. This ligated construct was then transformed into electrocompetent
DH5a cells (New England Biolabs) with time constants >5.1 at 1700mV. Cells
recovered for 45min at 30C and were then plated onto 7 different LB+kanamycin
(25 micrograms/mL) oversized plates (one for each final strain).

The resulting colonies were scraped and 800M cells from each plate were
inoculated into 200mL LB+kanamycin and grown at 30C overnight.

The plasmids were then isolated separately via Qiagen Maxiprep and cleaning
with a Zymo Clean and Concentrate kit (eluted in 12µL water).

The 7 different pLib-RiboJ::sfGFP plasmids were transformed into the 7
different strains (see Materials above). After 24h growth at 30C, these cells
were again scraped and 100M cells of each type were pooled together. 

To genome-integrate the plasmids at the galK locus, we previously made 7
strains (HG105, HG104, etc) with mCherry/chlor landing pads for AraC-inducible
Cre recombination at the galK locus. We followed a 5-day protocol to
genome-integrate the pLib libraries, which is performed as follows:

- **Day 1**: Grow glycerol stocks of each landing pad strains with the
  integration plasmid overnight in 200mL LB + kan (25ug/ml) at 30C.

- **Day 2**: Reinoculate culture with LB+0.2% Arabinose. Inoculate 200-400M
  cells from overnight cultures into 250mL LB + 0.2% arabinose and 25 (µg/mL)
  kan at 30C for 24 hours.

- **Day 3**: Pre-warm 4 LB + kan plates at 42C. Inoculate 800M cells of induced
  overnights into 80mL LB + kan at 42C. Grow cells to OD 0.3-0.7 (takes about
  1.5 hours). Plate 200M cells from the log-phase culture onto two plates, and
  grow for 16 hours at 42C. This heat-cures the pLib plasmids.

- **Day 4**: Check plates to ensure that there are enough colonies for 10X
  coverage.Scrape plates into 6mL LB. Inoculated 100M cells from EACH of the 7
  STRAINS into 200mL LB + kan. Grew overnight at 37C.

- **Day 5**: Make glycerol stocks of integrated libraries.

The pooled glycerol stocks were then inoculated overnight in LB with kan, to
saturation (this was done with two flasks in parallel). The next day, cells
were diluted 1:1000 into M9 media without antibiotic and with 0.5% glucose.
They were grown for 8 hours and then each culture was divided into three 50ml
falcons and three 15ml falcons and centrifuged at 4000xg for 10min. Pellets
were stored away. 50mL pellets were used for RNA extractions and 15mL pellets
were used to extract gDNA.

RNA was extracted from 50 mL library pellets using a Qiagen RNEasy Midi kit
(#75142) and 45 ug of each extract was concentrated using a Qiagen Minelute
Cleanup Kit (#74204). Barcoded cDNA was generated from 25 ug of each
concentrated RNA extract using Thermo Fisher SuperScript IV (#18090010) primed
with GU101. The manufacturer’s protocol was followed aside from extending the
reaction time to 1 hour at 52 °C. The cDNA reaction was cleaned using a Zymo
Research DNA Clean and Concentrator kit (#D40140) before amplification.

gDNA was extracted from 5 mL cell library pellets using a Qiagen Gentra
Puregene kit (#158567). Barcoded DNA was amplified from 1 ug of gDNA via PCR
for 14 cycles using primers GU59 and GU60. The reaction was subsequently
cleaned using a Zymo Research DNA Clean and Concentrator kit.

Sequencing adapters and indices were added as follows:

**PCR 1 - gDNA**

Perform qPCR to identify number of cycles necessary qPCR reaction (50 uL): 2x
Luna qPCR Mix (NEB): 25 µL H20: (23 µL - volume DNA) gDNA: 1 ug of DNA total
Primers (GU 59 + GU 60)@ 5 uM: 2 uL

16 cycles was mid-exponential

Repeat PCR with high fidelity polymerase (Q5). Performed 4 replicates for each
gDNA sample to have enough amplicon going forward. Amplify for two fewer cycles
than the number of cycles needed to reach mid-exponential phase as determined
by qPCR. Had 16 total tubes for this, 8 for each biological replicate, and 4
for each technical replicate. I wanted to start technical replicates in the
FIRST ROUND of PCR because this is where the bias will be introduced, if any.

Pooled each of the replicates together and performed PCR cleanup with the Zymo
Clean and Concentrator kit

**PCR 2 - gDNA**

Performed qPCR to identify number of cycles necessary qPCR reaction (50 uL): 2x
Luna qPCR Mix (NEB): 25 µL H20: (23 µL - volume DNA) PCR 1 amplicon: 2 µL @ 1
ng/µL Primers: 5 uM (GU 70 + GU 65-68)

9 cycles was mid-exponential

Repeat PCR with high fidelity polymerase (Q5). Perform 4 replicates for each
sample to have enough amplicon going forward. Amplify for two fewer cycles than
the number of cycles needed to reach mid-exponential phase as determined by
qPCR. Pool each of the 4 replicates together and perform PCR cleanup Make sure
you have ~300 bp product on a gel.

**PCR 1 - cDNA**

Perform qPCR to identify number of cycles necessary qPCR reaction (50 uL): 2x
Luna qPCR Mix (NEB): 25 µL H20: 21 µL cDNA: 2 µL (160ng) Primers: 5 uM (GU 59 +
GU 102) 16 cycles was mid-exponential

Compare number of cycles needed to get to mid-log in sample vs no reverse
transcriptase control. Ideally you’ll have >5 cycle difference before hitting
mid-log phase. There was no amplification of a reverse transcriptase control
after 20 cycles.

Repeat PCR with high fidelity polymerase (Q5). Perform 4 replicates for each
cDNA sample to have enough amplicon going forward. Amplify for two fewer cycles
than the number of cycles needed to reach mid-exponential phase as determined
by qPCR. Had 16 total tubes for this, 8 for each biological replicate, and 4
for each technical replicate. I wanted to start technical replicates in the
FIRST ROUND of PCR because this is where the bias will be introduced, if any.

Pooled each of the 4 replicates together and perform PCR cleanup

**PCR 2 - cDNA**

Perform qPCR to identify number of cycles necessary qPCR reaction (50 uL): 2x
Luna qPCR Mix (NEB): 25 µL H20: (23 µL - volume DNA) PCR 1 amplicon: 2 µL @ 1
ng/µL Primers: 5 uM (GU 102 + GU 61-64)

9 cycles was mid-exponential

Repeat PCR with high fidelity polymerase (Q5). Perform 4 replicates for each
sample to have enough amplicon going forward. Amplify for two fewer cycles than
the number of cycles needed to reach mid-exponential phase as determined by
qPCR. Pool each of the 4 replicates together and perform PCR cleanup Make sure
you have ~300 bp product on a gel.

ALL PCRs were performed with the following settings: 98C for 45s, and then X
cycles (16 or 9) of: 98C for 20s, 60C for 30s, 72C for 30s, and close with 72C
for 60s

**QC samples for Sequencing preparation:** Measure DNA concentration using
NanoDrop

Dilute samples to 20 nM (roughly 4 ng/uL) and pool in equimolar amounts. Make
sure the total volume of the pooled samples is greater than 20 uL.

Measure concentration using Qubit (final concentration: 4.17nM)

When I pooled the samples, I pooled them based on ng/µL concentrations, and
essentially added 200ng of each of the 8 libraries to a single tube in water.

GU61-GU62 are technical replicates of the first cDNA biological replicate
GU63-GU64 are technical replicates of the seconds cDNA biological replicate
GU65-GU66 are technical replicates of the first gDNA biological replicate
GU67-GU68 are technical replicates of the second gDNA biological replicate

For the sequencing, we spiked in 6.8µL of 50uM GU60 (Read 1 Primer) into
Illumina Well 12 We spiked in 6.8µL of 50uM GU71 (Index 1 Primer) into Illumina
Well 13

151 cycles for Read1, 50 cycles for Read2