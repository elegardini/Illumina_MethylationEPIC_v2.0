# =============================================================================
# HPA Axis Candidate Genes - CpG Extraction from EPICv2 Array
# Script: 01_extract_cpg_from_genes.R
# Description: Extract all CpG probes targeting HPA axis-related genes
#              from the Illumina EPICv2 manifest, including promoter regions
#              (10kb upstream). Also prepares files for EWAS Catalog queries.
# Author: Elena Gardini
# Date: 2026-04-24
# =============================================================================

# =============================================================================
# 1. LOAD LIBRARIES
# =============================================================================

library(biomaRt)
library(IlluminaHumanMethylationEPICv2anno.20a1.hg38)

# =============================================================================
# 2. DEFINE HPA AXIS GENES (36 genes from 7 functional groups)
# =============================================================================

all_genes <- c("CRH", "MC2R", "CRHR1", "CRHR2", "POMC",
               "NR3C1", "NR3C2", "FKBP5", "HSP90AA1",
               "CYP11B1", "CYP17A1", "CYP21A2", "HSD11B1", "HSD11B2",
               "SLC6A4", "HTR1A", "HTR2A", "COMT", "DRD2", "DRD3", 
               "GABRA1", "GABRB2", "GRIN2A", "GRM3",
               "BDNF", "NTRK2", "CREB1", "EGR1",
               "CRP", "IL6", "TNF",
               "CLOCK", "ARNTL", "CRY1", "CRY2", "PER1", "PER2")

cat("Total genes to analyze:", length(all_genes), "\n")

# =============================================================================
# 3. RETRIEVE GENOMIC COORDINATES FROM ENSEMBL (hg38)
# =============================================================================

# Connect to ENSEMBL database
ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")

# Download base coordinates
gene_coords <- getBM(
  attributes = c("hgnc_symbol", "chromosome_name", 
                 "start_position", "end_position"),
  filters = "hgnc_symbol",
  values = all_genes,
  mart = ensembl
)

# Keep only standard chromosomes (1-22, X, Y)
gene_coords <- gene_coords[gene_coords$chromosome_name %in% c(1:22, "X", "Y"), ]

# =============================================================================
# 4. EXTEND COORDINATES TO INCLUDE PROMOTER REGION (10KB UPSTREAM)
# =============================================================================

gene_coords_extended <- gene_coords
gene_coords_extended$start_extended <- gene_coords_extended$start_position - 10000
gene_coords_extended$start_extended[gene_coords_extended$start_extended < 0] <- 0
gene_coords_extended$end_extended <- gene_coords_extended$end_position + 10000

# Add "chr" prefix to match EPICv2 annotation format
gene_coords_extended$chr <- paste0("chr", gene_coords_extended$chromosome_name)

# Select columns for final coordinates table
gene_coords_final <- gene_coords_extended[, c("hgnc_symbol", "chr", "start_position", "end_position", 
                                              "start_extended", "end_extended")]

cat("\nFirst 5 genes with extended coordinates:\n")
print(head(gene_coords_final, 5))

# =============================================================================
# 5. EXTRACT CpG PROBES FROM EPICv2 ANNOTATION
# =============================================================================

# Load EPICv2 annotation
ann_epicv2 <- getAnnotation(IlluminaHumanMethylationEPICv2anno.20a1.hg38)

# Create empty dataframe with same structure as ann_epicv2
all_extended_probes <- ann_epicv2[0, ]
all_extended_probes$gene_name <- character()
all_extended_probes$region_type <- character()

# Loop through each gene and extract CpGs in extended region
for(i in 1:nrow(gene_coords_final)) {
  
  gene <- gene_coords_final[i, ]
  
  probes_in_region <- ann_epicv2[ann_epicv2$chr == gene$chr & 
                                   ann_epicv2$pos >= gene$start_extended & 
                                   ann_epicv2$pos <= gene$end_extended, ]
  
  if(nrow(probes_in_region) > 0) {
    probes_in_region$gene_name <- gene$hgnc_symbol
    probes_in_region$region_type <- "gene_promoter_10kb"
    all_extended_probes <- rbind(all_extended_probes, probes_in_region)
    cat("[OK]", gene$hgnc_symbol, "->", nrow(probes_in_region), "probes\n")
  } else {
    cat("[WARNING]", gene$hgnc_symbol, "-> NO probes found\n")
  }
}

# Remove duplicate rows
all_extended_probes <- unique(all_extended_probes)

# =============================================================================
# 6. SUMMARY STATISTICS
# =============================================================================

cat("\n=== TOTAL EXTENDED PROBES ===\n")
cat("Total probes:", nrow(all_extended_probes), "\n")
cat("Unique genes:", length(unique(all_extended_probes$gene_name)), "\n")

# =============================================================================
# 7. SAVE OUTPUT FILES
# =============================================================================

# 7a. Full CpG list for EWAS Catalog (one CpG per line)
cpg_list <- unique(all_extended_probes$Name)
writeLines(cpg_list, con = "HPA_Extended_CpG_list_for_EWAS.txt")
cat("\n[SAVED] CpG list: HPA_Extended_CpG_list_for_EWAS.txt (", length(cpg_list), " CpGs)\n")

# 7b. Reference table with CpG and gene information
cpg_reference_table <- all_extended_probes[, c("Name", "gene_name", "chr", "pos", 
                                               "Relation_to_Island", "UCSC_RefGene_Group")]
cpg_reference_table <- unique(cpg_reference_table)
cpg_reference_table <- cpg_reference_table[order(cpg_reference_table$gene_name), ]

write.csv(cpg_reference_table, 
          file = "HPA_Extended_CpG_with_Genes.csv", 
          row.names = FALSE)

write.table(cpg_reference_table[, c("Name", "gene_name")], 
            file = "HPA_Extended_CpG_list_with_Genes.txt", 
            sep = "\t", row.names = FALSE, quote = FALSE)

cat("[SAVED] HPA_Extended_CpG_with_Genes.csv (complete table)\n")
cat("[SAVED] HPA_Extended_CpG_list_with_Genes.txt (CpG + gene only)\n")
cat("\nTotal unique CpGs:", nrow(cpg_reference_table), "\n")

# =============================================================================
# 8. NOTE: For EWAS Catalog queries, use the API or website
# =============================================================================

cat("\n=== NEXT STEPS ===\n")
cat("To query EWAS Catalog:\n")
cat("1. Go to https://www.ewascatalog.org/\n")
cat("2. Upload HPA_Extended_CpG_list_for_EWAS.txt for batch search\n")
cat("3. Or search individually by gene name (e.g., BDNF, NR3C1)\n")

# =============================================================================
# End of script
# =============================================================================