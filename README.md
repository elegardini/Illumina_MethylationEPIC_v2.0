## Data Source
- GEO Accession: GSE229715
- Platform: Illumina MethylationEPIC v2.0 (GPL33022)
- Samples: 16 brain tumors of various histological types
- Original study: Turakulov R. et al. (2023)

## Requirements
- R (version 4.0 or higher)
- Bioconductor packages:
  - minfi
  - limma
  - DMRcate
  - IlluminaHumanMethylationEPICv2manifest
  - IlluminaHumanMethylationEPICv2anno.20a1.hg38
  - GEOquery
  - enrichR

## Key Results Summary
| Analysis | Result |
|----------|--------|
| DMRs identified | 5,852 age-associated regions |
| DMRs located in genes | 80.4% |
| Top genes | FOXP1, ZFHX3, RAD51B, HDAC4, CAVIN2 |
| PC1 variance | 18.5% |
| PC2 age association | 21.2% (R²) |

Due to the limited sample size (n=16) and the diversity of tumor types (each sample representing a different histological subtype), it is not possible to draw definitive conclusions about age-specific epigenetic alterations. The observed patterns should be interpreted with caution. The primary aim of this work is educational and methodological.

**Preliminary observations:**
- Age-associated DMRs were identified in genes involved in transcriptional regulation and DNA repair (FOXP1, ZFHX3, RAD51B, HDAC4, CAVIN2)
- However, these associations may be confounded by tumor type, as different tumors occur at different ages
- The analysis pipeline is reproducible and complete, but biological interpretations remain hypothesis-generating

## License
MIT License

## Author
Elena Gardini

## Contact
elegardini@gmail.com
