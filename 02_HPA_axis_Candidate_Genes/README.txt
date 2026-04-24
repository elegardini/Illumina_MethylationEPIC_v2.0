# HPA Axis Candidate Genes - CpG Extraction from EPICv2 Array

## Objective
Extract all EPICv2 CpG probes targeting HPA axis-related genes, including promoter regions (10kb upstream), to enable targeted analysis of stress-response and psychiatric disorder-related epigenetic marks.

## Gene List (36 genes)

| Category | Genes |
|----------|-------|
| Activation feedback | CRH, MC2R, CRHR1, CRHR2, POMC |
| Glucocorticoid signaling | NR3C1, NR3C2, FKBP5, HSP90AA1 |
| Cortisol metabolism | CYP11B1, CYP17A1, CYP21A2, HSD11B1, HSD11B2 |
| Brain circuit | SLC6A4, HTR1A, HTR2A, COMT, DRD2, DRD3, GABRA1, GABRB2, GRIN2A, GRM3 |
| Neuroplasticity | BDNF, NTRK2, CREB1, EGR1 |
| Immune-HPA | CRP, IL6, TNF |
| Circadian rhythm | CLOCK, ARNTL, CRY1, CRY2, PER1, PER2 |

## Methods

1. **Retrieve genomic coordinates** for each gene using `biomaRt` (hg38)
2. **Extend coordinates** to include 10kb promoter region upstream
3. **Extract CpG probes** from EPICv2 annotation (`IlluminaHumanMethylationEPICv2anno.20a1.hg38`)
4. **Annotate** with CpG island relation and gene information

## Results

| Metric | Value |
|--------|-------|
| Total probes extracted | **1,896** |
| Unique genes covered | 36 |
| Probe range per gene | 12 (CRP) to 106 (GRIN2A) |

## Output Files

| File | Description |
|------|-------------|
| `HPA_ALL_GENES_EXTENDED_PROMOTER_10kb.csv` | Complete probe annotation (all columns) |
| `HPA_Extended_CpG_with_Genes.csv` | Probe-gene mapping with CpG island info |
| `HPA_Extended_CpG_list_for_EWAS.txt` | CpG list for EWAS Catalog queries |

## Script

Run `scripts/01_extract_cpg_from_genes.R` to reproduce the analysis.

## Requirements

- R (version 4.5.2)
- Bioconductor packages:
  - `biomaRt`
  - `IlluminaHumanMethylationEPICv2anno.20a1.hg38`

## Author
Elena Gardini