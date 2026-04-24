# Illumina MethylationEPIC v2.0 Array Analysis

This repository contains two independent projects analyzing DNA methylation data.

## Projects

### 1. Brain Tumor Aging Analysis
Age-associated methylation changes in 16 brain tumor samples (GSE229715).
➡️ [View project](01_Brain_Tumor_Aging/README.md)

### 2. HPA Axis Candidate Genes - CpG Extraction

This section demonstrates how to extract all CpG probes targeting HPA axis-related genes 
from the Illumina EPICv2 manifest, including promoter regions (10kb upstream).

**Genes analyzed**: CRH, MC2R, CRHR1, CRHR2, POMC, NR3C1, NR3C2, FKBP5, HSP90AA1,
CYP11B1, CYP17A1, CYP21A2, HSD11B1, HSD11B2, SLC6A4, HTR1A, HTR2A, COMT, DRD2, 
DRD3, GABRA1, GABRB2, GRIN2A, GRM3, BDNF, NTRK2, CREB1, EGR1, CRP, IL6, TNF, 
CLOCK, ARNTL, CRY1, CRY2, PER1, PER2.

**Output**: 1,896 unique CpG probes with genomic coordinates, CpG island relations, 
and gene annotations.

➡️ [View HPA analysis details](02_HPA_Axis_Candidate_Genes/README.md)

## Author
Elena Gardini

## License
MIT
