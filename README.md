# Illumina MethylationEPIC v2.0 Array Analysis

This repository contains two independent projects analyzing DNA methylation data from the Illumina EPICv2 (950K) array.

## Projects

### 1. Brain Tumor Aging Analysis
Age-associated methylation changes in 16 brain tumor samples (GSE229715).

**Key findings**: Age-associated DMRs in genes involved in transcriptional regulation and DNA repair (FOXP1, ZFHX3, RAD51B, HDAC4, CAVIN2).

➡️ [View project details](01_Brain_Tumor_Aging/README.md)

### 2. HPA Axis Related Candidate Genes - CpG Extraction

Extraction of all CpG probes targeting HPA axis / HPA axis-related genes from the Illumina EPICv2 manifest, including promoter regions (10kb upstream).

**Genes analyzed**: 36 HPA axis-related genes including NR3C1, BDNF, FKBP5, SLC6A4, CRH, and others.

**Output**: 1,896 unique CpG probes with genomic coordinates, CpG island relations, and gene annotations.

➡️ [View project details](02_HPA_axis_Candidate_Genes/README.md)

## Author
Elena Gardini

## License
MIT
