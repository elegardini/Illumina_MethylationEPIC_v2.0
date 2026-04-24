
#Identify and download dataset: GSE229715 - Brain Tumor Classification- 16 samples. Description: Study on brain tumor classification comparing EPIC v1 (850K) and v2 (950K) arrays.Includes both EPIC v1 and v2 data, allowing direct comparison between platforms
#Download dataset and metadata: In R. 
library(GEOquery)
metadata <- "C:/Users/elega/OneDrive/Desktop/EPIC Project/metadata"
gse_id <- "GSE229715"
gse <- getGEO(gse_id, destdir = metadata)
#File Epic v2 (GPL33022)
gse_v2 <- getGEO(filename = "metadata/GSE229715-GPL33022_series_matrix.txt.gz")
#Samples information (pheno_data)
pheno_data <- pData(phenoData(gse_v2))
print(pheno_data[1:3, "characteristics_ch1"])
char_cols <- grep("characteristics", colnames(pheno_data), value = TRUE)
print(char_cols)
#all characteristics for sample 1
cat("\nAll info first sample:\n")
for(col in char_cols) {
  cat(col, ":", pheno_data[1, col], "\n")
}
# Create a sample sheet with chosen information (sample name, disease, sex,age,material, tumor site)
sample_info <- data.frame(
  Sample_Name = pheno_data$geo_accession,  # or pheno_data$title
  Disease = NA,
  Sex = NA,
  Age = NA,
  Material = NA,
  Tumor_Site = NA
)

# Extract information from all characteristics columns
char_cols <- grep("characteristics", colnames(pheno_data), value = TRUE)

for(i in 1:nrow(pheno_data)) {
  for(col in char_cols) {
    value <- pheno_data[i, col]
    if(!is.na(value)) {
      if(grepl("disease", value, ignore.case = TRUE)) {
        sample_info$Disease[i] <- gsub(".*disease state:\\s*", "", value)
      }
      if(grepl("Sex", value, ignore.case = TRUE)) {
        sample_info$Sex[i] <- gsub(".*Sex:\\s*", "", value)
      }
      if(grepl("age", value, ignore.case = TRUE)) {
        sample_info$Age[i] <- gsub(".*age:\\s*", "", value)
      }
      if(grepl("material", value, ignore.case = TRUE)) {
        sample_info$Material[i] <- gsub(".*material_type:\\s*", "", value)
      }
      if(grepl("tumor_site", value, ignore.case = TRUE)) {
        sample_info$Tumor_Site[i] <- gsub(".*tumor_site:\\s*", "", value)
      }
    }
  }
}

# View the result
cat("Sample sheet created:\n")
print(sample_info)

# Save file
write.csv(sample_info, "metadata/sample_sheet_epicv2.csv", row.names = FALSE)

#extracting row methylation data
#idata file: we extract the intensity signals red and green (it s entire number like, 11145 or 450), then the 2 intensities red and green are matched toghether, there is then a preprocessing normalization (per correggere i batch effect, it's a calibration to make 
things comparable, for example the Noob (Normal-exponential out of band) correct the background noise, after that we transform the signals in methylation values Beta or M-values (which for each site is the rate methylated/methylated+non-methylated).

#Download idat files
#extract the GSM IDs 
sample_info <- read.csv("Metadata_sample_sheet_epicv2_GPL33022.csv")
gsm_ids <- sample_info$Sample_Name
print(gsm_ids)
# download all IDAT files
library(GEOquery)
getwd()
idat_dir<-"C:/Users/elega/OneDrive/Desktop/EPIC Project/data_idat"
for(gsm in gsm_ids) {
  cat("Downloading:", gsm, "...\n")
  try(getGEOSuppFiles(gsm, baseDir = idat_dir))
}
#Load IDAT files with minfi (EPIC v2)
library(minfi)
library(IlluminaHumanMethylationEPICv2manifest)
library(IlluminaHumanMethylationEPICv2anno.20a1.hg38)
base_dir <- "C:/Users/elega/OneDrive/Desktop/EPIC Project"
idat_dir <- file.path(base_dir, "data_idat")
metadata_dir <- file.path(base_dir, "metadata")
# Read metadata sheet
sample_info <- read.csv(file.path(metadata_dir, "Metadata_sample_sheet_epicv2_GPL33022.csv"))
print(head(sample_info))
# Get all GSM folders
gsm_folders <- list.files(idat_dir, pattern = "GSM", full.names = TRUE)
cat("Found", length(gsm_folders), "GSM folders\n")
# Create targets dataframe (map linking samples to their IDAT files and clinical data)
targets <- data.frame(
  Sample_Name = sample_info$Sample_Name,
  Basename = file.path(idat_dir, sample_info$Sample_Name, 
                       sample_info$Sample_Name),  # Path to IDAT files (without _Grn/_Red)
  Disease = sample_info$Disease,
  Sex = sample_info$Sex,
  Age = sample_info$Age,
  Material = sample_info$Material,
  Tumor_Site = sample_info$Tumor_Site,
  stringsAsFactors = FALSE
)
#check
cat("Targets map created successfully!\n")
cat("First 3 rows of the targets map:\n")
print(targets[1:3, 1:4])  # Show first 3 rows and first 4 columns
#extract methylation data
# Decompress all .gz files before loading
idat_dir <- "C:/Users/elega/OneDrive/Desktop/EPIC Project/data_idat"
# Find all .gz files
gz_files <- list.files(idat_dir, pattern = "\\.gz$", recursive = TRUE, full.names = TRUE)
cat("Found", length(gz_files), "compressed files\n")

# Decompress each file
for(gz_file in gz_files) {
  cat("Decompressing:", basename(gz_file), "\n")
  R.utils::gunzip(gz_file, remove = FALSE)  # remove = FALSE keeps the original .gz
}

# Rename all files by removing the extra part
idat_dir <- "C:/Users/elega/OneDrive/Desktop/EPIC Project/data_idat"

# Find all .idat files
idat_files <- list.files(idat_dir, pattern = "\\.idat$", recursive = TRUE, full.names = TRUE)

for(file in idat_files) {
  # Extract the file name
  old_name <- basename(file)
  
  # Create the new name (only GSM_Grn.idat or GSM_Red.idat)
  # Example: "GSM7173814_207219640008_R01C01_Grn.idat" -> "GSM7173814_Grn.idat"
  parts <- strsplit(old_name, "_")[[1]]
  new_name <- paste(parts[1], parts[length(parts)], sep = "_")
  
  # Rename the file
  new_path <- file.path(dirname(file), new_name)
  file.rename(file, new_path)
  cat("Renamed:", old_name, "->", new_name, "\n")
}

#load the row methylation data using the target map-pre-normalization-separate green and red signals
rgSet <- read.metharray.exp(targets = targets, force = TRUE)
#check
rgSet
cat("=== RGChannelSet SUMMARY ===\n")
cat("Number of samples:", ncol(rgSet), "\n")
cat("Number of probes:", nrow(rgSet), "\n")
cat("Class:", class(rgSet), "\n")
cat("\nSample names:\n")
print(sampleNames(rgSet))
cat("\nFirst rows of green channel (first 2 samples):\n")
print(head(getGreen(rgSet)[,1:2]))
#save
save(rgSet, targets, file = "C:/Users/elega/OneDrive/Desktop/EPIC Project/data_processed/EPICv2_rgSet_complete.RData")

#Normalization and Beta values. Normalizazion using Noob. What it does is Background correction (removes noise),Dye-bias correction (balances red/green channels),Calculates methylation values (combines red+green into one value per CpG), these methylation values are still intensities
Beta values are obtained using the getBeta function
# Normalization with Noob 2026-03-04 Load raw data
load("C:/Users/elega/OneDrive/Desktop/EPIC Project/data_processed/EPICv2_rgSet_complete.RData")
library(minfi)
#check
cat("Data loaded!\n")
cat("Samples:", ncol(rgSet), "\n")
cat("Probes:", nrow(rgSet), "\n")
#running Noob normalization
cat("Running Noob normalization...\n")
grSet <- preprocessNoob(rgSet)
#check
cat("Normalization complete!\n")
cat("Class of grSet:", class(grSet), "\n")
cat("Number of samples:", ncol(grSet), "\n")
cat("Number of probes:", nrow(grSet), "\n")
cat("\nSample names:\n")
print(sampleNames(grSet)[1:3])
#extract beta values 
beta <- getBeta(grSet)
#check 
cat("Beta matrix dimensions:", dim(beta), "\n")
cat("Beta value range:", range(beta, na.rm = TRUE), "\n")
cat("\nFirst 5 probes for first 3 samples:\n")
print(beta[1:5, 1:3])
cat("\nSummary of first sample (GSM7173814):\n")
summary(beta[,1])
#save beta values file
save(grSet, beta, targets, 
     file = "C:/Users/elega/OneDrive/Desktop/EPIC Project/data_processed/EPICv2_beta_values.RData")

#Density plot (quality control for distribution-biomodal=normal distribution for metylation data), PCA (principal component analysis-general) to find the CpGs patterns explaining most of the variance, we identify the top 2 PC, PC1(18.5%) and PC2 (14.8%),
#then PCA colored by clinical variables: where mixed colours along PC1 and PC2 = the variable doesn't explain the variance, where colors are separated= variable somhow explained. Age was the only one who looked like explain part of the varaince, so only for this variable we quantified the age effect
#and found that age explained 21.2% of PC2. 
#Quality control
# Load beta values
load("C:/Users/elega/OneDrive/Desktop/EPIC Project/data_processed/EPICv2_beta_values.RData")
#Plot first sample 
plot(density(beta[,1], na.rm = TRUE), 
     col = 1, 
     main = "Beta values distribution: first sample",
     xlab = "Beta", 
     ylab = "Density",
     ylim = c(0, 5))
# Add samples 2 through 16
for(i in 2:ncol(beta)) {
  lines(density(beta[,i], na.rm = TRUE), col = i)
}
# Legend with first 5 samples
legend("topright", legend = colnames(beta)[1:5], col = 1:5, lty = 1, cex = 0.6)
# Run PCA (principal components analysis) on beta values
pca <- prcomp(t(na.omit(beta)), center = TRUE, scale. = FALSE)
plot(pca$x[,1], pca$x[,2], 
     pch = 19,
     xlab = paste0("PC1 (", round(summary(pca)$importance[2,1]*100, 1), "%)"),
     ylab = paste0("PC2 (", round(summary(pca)$importance[2,2]*100, 1), "%)"),
     main = "PCA of 16 EPICv2 samples (no colors)")
# Add sample numbers
text(pca$x[,1], pca$x[,2], labels = 1:16, pos = 3, cex = 0.7)

#PCA matched with sample info (is the variation explained by phenotypes, if yes, separated colors, if not mixed colors)
sample_info <- read.csv("C:/Users/elega/OneDrive/Desktop/EPIC Project/metadata/Metadata_sample_sheet_epicv2_GPL33022.csv")

#create colors for Disease variable
disease_colors <- as.numeric(factor(sample_info$Disease))
# Plot
plot(pca$x[,1], pca$x[,2], 
     col = disease_colors, pch = 19, cex = 1.2,
     xlab = paste0("PC1 (", round(summary(pca)$importance[2,1]*100, 1), "%)"),
     ylab = paste0("PC2 (", round(summary(pca)$importance[2,2]*100, 1), "%)"),
     main = "PCA colored by Disease")

legend("topright", legend = unique(sample_info$Disease), 
       col = 1:length(unique(sample_info$Disease)), pch = 19)
# Second variable: sex. Create colors for Sex
sex_colors <- as.numeric(factor(sample_info$Sex))

# Plot
plot(pca$x[,1], pca$x[,2], 
     col = sex_colors, pch = 19, cex = 1.2,
     xlab = paste0("PC1 (", round(summary(pca)$importance[2,1]*100, 1), "%)"),
     ylab = paste0("PC2 (", round(summary(pca)$importance[2,2]*100, 1), "%)"),
     main = "PCA colored by Sex")

legend("topright", legend = unique(sample_info$Sex), 
       col = 1:length(unique(sample_info$Sex)), pch = 19)

#Third variable: age 
cat("Age values:\n")
print(sample_info$Age)
#convert 5M in 5 months
sample_info$Age_clean <- sample_info$Age
sample_info$Age_clean[sample_info$Age == "5M"] <- 0.42
sample_info$Age_numeric <- as.numeric(sample_info$Age_clean)
cat("Age values (original):", sample_info$Age, "\n")
cat("Age values (numeric):", sample_info$Age_numeric, "\n")

# Create color gradient by Age (5 groups)
age_colors <- as.numeric(cut(sample_info$Age_numeric, breaks = 5))

# Plot
plot(pca$x[,1], pca$x[,2], 
     col = age_colors, pch = 19, cex = 1.2,
     xlab = paste0("PC1 (", round(summary(pca)$importance[2,1]*100, 1), "%)"),
     ylab = paste0("PC2 (", round(summary(pca)$importance[2,2]*100, 1), "%)"),
     main = "PCA colored by Age (cleaned)")

legend("topright", 
       legend = levels(cut(sample_info$Age_numeric, breaks = 5)), 
       col = 1:5, pch = 19, cex = 0.7)
# QUANTIFY HOW MUCH VARIANCE IS EXPLAINED BY AGE

# 1. Linear model: PC1 ~ Age
lm_pc1 <- lm(pca$x[,1] ~ sample_info$Age_numeric)
summary_pc1 <- summary(lm_pc1)
cat("R² for PC1 ~ Age:", round(summary_pc1$r.squared * 100, 1), "%\n")
# 2. Linear model: PC2 ~ Age (result: it explains 21.2% of the variance of CP2)
lm_pc2 <- lm(pca$x[,2] ~ sample_info$Age_numeric)
summary_pc2 <- summary(lm_pc2)
cat("R² for PC2 ~ Age:", round(summary_pc2$r.squared * 100, 1), "%\n")
#fourth variable: material 
# Check unique values in Material
cat("Unique Material values:\n")
print(unique(sample_info$Material))
#fifth variable: tumor site
# Check unique Tumor_Site values
cat("Unique Tumor_Site values:\n")
print(unique(sample_info$Tumor_Site))

#PC1 Loadings Analysis (which CpG sites are driving the main variation (PC1 = 18.5%)), DMP analyses (not finished because no FDR significant),  enrichement ( GO Enrichment Analysis on age-associated CpGs (top 10,000 from DMP model), revealing enrichment in development genes (HOXA9, HOXB8, MYB, etc.)
#PC1 loadings (2026-03-10)
# Load beta values
library(minfi)
load("C:/Users/elega/OneDrive/Desktop/EPIC Project/data_processed/EPICv2_beta_values.RData")
#run PCA-analysis for PCA to recall in memory
pca <- prcomp(t(na.omit(beta)), center = TRUE, scale. = FALSE)
# Extract loadings for PC1 (weights for each CpG)
loadings_pc1 <- pca$rotation[,1]  
# Check
cat("PC1 loadings:\n")
cat("Number of CpGs:", length(loadings_pc1), "\n")
cat("Range:", range(loadings_pc1), "\n")
# Top 10 CpGs with largest POSITIVE loadings
top_positive <- sort(loadings_pc1, decreasing = TRUE)[1:10]
cat("\nTop 10 positive loadings (drive high methylation in high PC1 samples):\n")
print(top_positive)
# Top 10 CpGs with largest NEGATIVE loadings  
top_negative <- sort(loadings_pc1, decreasing = FALSE)[1:10]
cat("\nTop 10 negative loadings (drive low methylation in high PC1 samples):\n")
print(top_negative)
# Load annotation package- to check in which gene/regions are located the top 20 cpgs and their function if known
library(IlluminaHumanMethylationEPICv2anno.20a1.hg38)
# Get annotation for all CpGs
annotation <- getAnnotation(IlluminaHumanMethylationEPICv2anno.20a1.hg38)
# Annotate top positive CpGs
cat("\n=== TOP POSITIVE CpGs (high methylation in high PC1 samples) ===\n")
pos_anno <- annotation[names(top_positive), 
                       c("chr", "pos", "UCSC_RefGene_Name", "Regulatory_Feature_Group")]
print(pos_anno)
# Annotate top negative CpGs
cat("\n=== TOP NEGATIVE CpGs (low methylation in high PC1 samples) ===\n")
neg_anno <- annotation[names(top_negative), 
                       c("chr", "pos", "UCSC_RefGene_Name", "Regulatory_Feature_Group")]
print(neg_anno)
#filtering and DMP analysis
#filtering= remove probes with signal too close to noise, SNP affected probes (Probes near genetic variants (SNPs) can give spurious methylation signals,based on distance to SNP (e.g., within 2bp) and minor allele frequency),probe sequences that bind to multiple genomic locations (according to listo of Chen 2013 and Pidsley 2016), sex chormosome if not studying sex differences
# Load libraries
library(DMRcate)
library(IlluminaHumanMethylationEPICv2anno.20a1.hg38)
#filtering, including sex chromosomes
#Run filtering (this removes SNPs, cross-hybridizing probes, AND sex chromosomes)
beta_filtered <- rmSNPandCH(beta, 
                            dist = 2,           # Remove SNPs within 2bp
                            mafcut = 0.05,      # Minor allele frequency cutoff
                            rmcrosshyb = TRUE,  # Remove cross-hybridizing probes
                            rmXY = TRUE)        # Remove X and Y chromosomes

# Check new dimensions
cat("Filtered beta dimensions:", dim(beta_filtered), "\n")
cat("Probes removed:", nrow(beta) - nrow(beta_filtered), "\n")
cat("Percentage kept:", round(nrow(beta_filtered)/nrow(beta)*100, 1), "%\n")
# Save filtered beta values with clear filename
save(beta_filtered, file = "C:/Users/elega/OneDrive/Desktop/EPIC Project/data_processed/EPICv2_beta_filtered_sexchr_removed.RData")
# Confirm we can load it back
load("C:/Users/elega/OneDrive/Desktop/EPIC Project/data_processed/EPICv2_beta_filtered_sexchr_removed.RData")
cat("Reloaded beta_filtered dimensions:", dim(beta_filtered), "\n")
#DMP analysis- only age as the most intersting and there is too much variation or 0 variation in other variables
## Load metadata
sample_info <- read.csv("C:/Users/elega/OneDrive/Desktop/EPIC Project/metadata/Metadata_sample_sheet_epicv2_GPL33022.csv")
# Clean Age (we did this before)
sample_info$Age_clean <- sample_info$Age
sample_info$Age_clean[sample_info$Age == "5M"] <- 0.42
sample_info$Age_numeric <- as.numeric(sample_info$Age_clean)
# Verify samples are in same order in the metadata file and in the beta filtered file
cat("Samples match?", all(colnames(beta_filtered) == sample_info$Sample_Name), "\n")
# Load limma for DMP analyses
library(limma)
# Design matrix: ~ Age (continuous)
design_age <- model.matrix(~ Age_numeric, data = sample_info)
# Fit model
fit <- lmFit(beta_filtered, design_age)
fit <- eBayes(fit)
# Results
age_dmp <- topTable(fit, coef = "Age_numeric", number = Inf, adjust.method = "BH")

# check
cat("CpGs associated with age (FDR < 0.05):", sum(age_dmp$adj.P.Val < 0.05), "\n")
cat("\nTop 10 age-associated CpGs:\n")
print(head(age_dmp, 10))
#Pathways analyses: see if with of string of cpgs we can detect a difference: with 16 samples, the FDR correction method above showed there not enough statistical power to detect individual CpG effects at a genome-wide significant level.
#so we do this with go enrichment 
# BiocManager::install("missMethyl")

library(missMethyl)
library(IlluminaHumanMethylationEPICv2anno.20a1.hg38)
# Get top 10,000 CpGs from age analysis
sig_cpgs <- rownames(age_dmp)[1:10000]

# Get all CpGs that were tested
all_cpgs <- rownames(beta_filtered)

cat("Number of significant CpGs (top 10,000):", length(sig_cpgs), "\n")
cat("Number of background CpGs:", length(all_cpgs), "\n")

# Run GO enrichment analysis= check if correspond to a certain biological function 
install.packages("enrichR")
library(enrichR)
# Get gene names for your top CpGs
library(IlluminaHumanMethylationEPICv2anno.20a1.hg38)
annotation <- getAnnotation(IlluminaHumanMethylationEPICv2anno.20a1.hg38)
# Map CpGs to genes
sig_genes <- unique(unlist(strsplit(annotation[sig_cpgs[1:1000], "UCSC_RefGene_Name"], ";")))
sig_genes <- sig_genes[sig_genes != ""]
#Run enrichment with enrichR
dbs <- c("GO_Biological_Process_2023", "GO_Molecular_Function_2023", "GO_Cellular_Component_2023")
enriched <- enrichr(sig_genes, dbs)
# View results
print(head(enriched$GO_Biological_Process_2023))
# Save just topGO 
save(topGO, file = "C:/Users/elega/OneDrive/Desktop/EPIC Project/results/GO_topTerms.RData")

# Also save the genes list
interesting_genes <- c("HOXA9", "HOXB8", "HOXA5", "MYB", "GREM1", "WT1", "NRG1", "ASB2", "EOMES", "TBX20")
save(interesting_genes, file = "C:/Users/elega/OneDrive/Desktop/EPIC Project/results/age_associated_genes.RData")

#DMR analysis with DMRcate- conclusion: Observation of age-associated DNA methylation changes in a set of brain tumor samples, with enrichment in genes involved in transcriptional regulation and DNA repair. 
#Among the most frequently affected genes were FOXP1, ZFHX3, RAD51B, HDAC4 (predominantly hypermethylated with age) and CAVIN2 (hypomethylated). 
#These genes have been previously associated with tumor suppression and cancer progression. Data suggest that aging is associated with epigenetic remodeling in key genes involved in tumor suppression and DNA repair. 
#This could potentially translate into enhanced tumor aggressiveness and diminished treatment response in older patients.

#DMR analysis with DMRcate
# Load filtered beta values
load("C:/Users/elega/OneDrive/Desktop/EPIC Project/data_processed/EPICv2_beta_filtered_sexchr_removed.RData")
# Load metadata
sample_info <- read.csv("C:/Users/elega/OneDrive/Desktop/EPIC Project/metadata/Metadata_sample_sheet_epicv2_GPL33022.csv")
# Clean Age
sample_info$Age_clean <- sample_info$Age
sample_info$Age_clean[sample_info$Age == "5M"] <- 0.42
sample_info$Age_numeric <- as.numeric(sample_info$Age_clean)
# Convert beta to M-values (DMRcate works better with M-values-log2 of M-values)
library(DMRcate)
M_vals <- logit2(beta_filtered)
# Create design matrix 
design_age <- model.matrix(~ Age_numeric, data = sample_info)
# Create annotation object for DMRcate- statistics for each CpG site variation in function of age
myannotation <- cpg.annotate("array", 
                             object = M_vals, 
                             what = "M",
                             arraytype = "EPICv2", 
                             analysis.type = "differential",
                             design = design_age,
                             coef = "Age_numeric",
                             fdr = 0.05)

# Run DMR detection. DMR associated with age
dmr_results <- dmrcate(myannotation, 
                       lambda = 1000,      # 1kb smoothing window
                       C = 2)               # Scaling factor
#Error in dmrcate(myannotation, lambda = 1000, C = 2) : The FDR specified in cpg.annotate() returned no significant CpGs, hence there are no DMRs.Try specifying avalue of 'pcutoff' in dmrcate() and/or increasing 'fdr' in cpg.annotate(). 
##no significan cpg with FDR, so we use manual cuoff=0.05= pvalue not corrected
# Run DMR detection with manual p-value cutoff
dmr_results <- dmrcate(myannotation, 
                       lambda = 1000,      
                       C = 2,
                       pcutoff = 0.05)      # Usa p-value grezzo, non FDR
# Extract results
dmr_ranges <- extractRanges(dmr_results, genome = "hg38")
# Check how many DMRs found
cat("Number of DMRs found:", length(dmr_ranges), "\n")
# Create and view DMR table
dmr_table <- data.frame(
  chr = seqnames(dmr_ranges),
  start = start(dmr_ranges),
  end = end(dmr_ranges),
  width = width(dmr_ranges),
  nCpGs = dmr_ranges$no.cpgs,
  minFDR = dmr_ranges$min_smoothed_fdr,
  Stouffer = dmr_ranges$Stouffer,
  maxdiff = dmr_ranges$maxdiff,
  meandiff = dmr_ranges$meandiff,
  overlapping_genes = dmr_ranges$overlapping.genes
)

# Look at first few
cat("\nTop 10 DMRs:\n")
print(head(dmr_table, 10))

# Summary statistics
cat("\n=== DMR SUMMARY ===\n")
cat("Median DMR width:", median(dmr_table$width), "bp\n")
cat("Median number of CpGs per DMR:", median(dmr_table$nCpGs), "\n")
cat("DMRs with gene annotation:", sum(dmr_table$overlapping_genes != ""), "\n")

# Check what's in overlapping_genes
cat("First few values:\n")
print(head(dmr_table$overlapping_genes))

# Count DMRs with gene annotation (handling NA)
dmr_table$overlapping_genes[is.na(dmr_table$overlapping_genes)] <- ""
cat("\nDMRs with gene annotation:", sum(dmr_table$overlapping_genes != ""), "\n")

# Percentage
cat("Percentage with genes:", round(sum(dmr_table$overlapping_genes != "")/nrow(dmr_table)*100, 1), "%\n")

# Find genes with most DMRs
library(stringr)
all_genes <- unlist(str_split(dmr_table$overlapping_genes[dmr_table$overlapping_genes != ""], ", "))
gene_freq <- sort(table(all_genes), decreasing = TRUE)
cat("\nTop 10 genes with most DMRs:\n")
print(head(gene_freq, 10))
# DMR with meandiff- hypo or hyper methylated with age change 
print(dmr_table[1:10, c("overlapping_genes", "meandiff")])
# Save
write.csv(dmr_table, 
          file = "C:/Users/elega/OneDrive/Desktop/EPIC Project/results/DMR_age_associated.csv",
          row.names = FALSE)
