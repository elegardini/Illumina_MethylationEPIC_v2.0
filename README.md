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

These findings suggest that aging in brain tumors may be associated with epigenetic alterations, including hypermethylation of genes implicated in transcriptional regulation and DNA repair. In particular, hypermethylation of FOXP1, ZFHX3, RAD51B, and HDAC4 could reflect age-related epigenetic remodeling that may influence tumor biology. Additionally, the observed changes in CAVIN2 expression may be consistent with pathways related to tumor invasiveness.

## License
MIT License

## Author
Elena Gardini

## Contact
elegardini@gmail.com
