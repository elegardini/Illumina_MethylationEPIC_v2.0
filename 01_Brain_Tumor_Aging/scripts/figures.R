# ====================================================
# Generate figures for GitHub portfolio
# Script: 05_generate_figures.R
# ====================================================

#packages
install.packages("pheatmap")
# Load libraries
library(ggplot2)
library(pheatmap)
library(RColorBrewer)
library(stringr)
library(minfi)

# Set output directory
fig_dir <- "C:/Users/elega/OneDrive/Desktop/EPIC Project/GitHub_Repo/figures"
dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)

# ====================================================
# 1. Load data
# ====================================================
cat("Loading data...\n")
load("C:/Users/elega/OneDrive/Desktop/EPIC Project/data_processed/EPICv2_beta_values.RData")
load("C:/Users/elega/OneDrive/Desktop/EPIC Project/data_processed/EPICv2_beta_filtered_sexchr_removed.RData")
sample_info <- read.csv("C:/Users/elega/OneDrive/Desktop/EPIC Project/metadata/Metadata_sample_sheet_epicv2_GPL33022.csv")
dmr_table <- read.csv("C:/Users/elega/OneDrive/Desktop/EPIC Project/results/DMR_age_associated.csv")

# Clean Age
sample_info$Age_clean <- sample_info$Age
sample_info$Age_clean[sample_info$Age == "5M"] <- 0.42
sample_info$Age_numeric <- as.numeric(sample_info$Age_clean)

# ====================================================
# 2. Figure 1: Density plot
# ====================================================
cat("Generating Figure 1: Density plot...\n")
png(file.path(fig_dir, "01_density_plot.png"), width = 800, height = 600, res = 100)

plot(density(beta[,1], na.rm = TRUE), 
     col = 1, lwd = 2,
     main = "Beta values distribution across 16 samples",
     xlab = "Beta (methylation level)", 
     ylab = "Density",
     ylim = c(0, 5))

for(i in 2:ncol(beta)) {
  lines(density(beta[,i], na.rm = TRUE), col = i, lwd = 1.5)
}
legend("topright", legend = colnames(beta)[1:5], 
       col = 1:5, lty = 1, lwd = 2, cex = 0.8)

dev.off()

# ====================================================
# 3. Figure 2: PCA (no colors)
# ====================================================
cat("Generating Figure 2: PCA plot...\n")
pca <- prcomp(t(na.omit(beta)), center = TRUE, scale. = FALSE)

png(file.path(fig_dir, "02_pca_nocolors.png"), width = 800, height = 600, res = 100)

plot(pca$x[,1], pca$x[,2], 
     pch = 19, cex = 1.2,
     xlab = paste0("PC1 (", round(summary(pca)$importance[2,1]*100, 1), "%)"),
     ylab = paste0("PC2 (", round(summary(pca)$importance[2,2]*100, 1), "%)"),
     main = "PCA of 16 EPICv2 samples")
text(pca$x[,1], pca$x[,2], labels = 1:16, pos = 3, cex = 0.8)

dev.off()

# ====================================================
# 4. Figure 3: PCA colored by Age
# ====================================================
cat("Generating Figure 3: PCA colored by Age...\n")
png(file.path(fig_dir, "03_pca_age.png"), width = 800, height = 600, res = 100)

age_colors <- as.numeric(cut(sample_info$Age_numeric, breaks = 5))
plot(pca$x[,1], pca$x[,2], 
     col = age_colors, pch = 19, cex = 1.5,
     xlab = paste0("PC1 (", round(summary(pca)$importance[2,1]*100, 1), "%)"),
     ylab = paste0("PC2 (", round(summary(pca)$importance[2,2]*100, 1), "%)"),
     main = "PCA colored by Age (quintiles)")

legend("topright", 
       legend = levels(cut(sample_info$Age_numeric, breaks = 5)), 
       col = 1:5, pch = 19, cex = 0.8)

dev.off()

# ====================================================
# 5. Figure 4: Top genes with most DMRs (barplot)
# ====================================================
cat("Generating Figure 4: Top genes barplot...\n")

# Extract top genes from DMR table
dmr_table$overlapping_genes[is.na(dmr_table$overlapping_genes)] <- ""
all_genes <- unlist(str_split(dmr_table$overlapping_genes[dmr_table$overlapping_genes != ""], ", "))
gene_freq <- sort(table(all_genes), decreasing = TRUE)
top10 <- head(gene_freq, 10)

png(file.path(fig_dir, "04_top_genes_dmrs.png"), width = 900, height = 600, res = 100)

par(mar = c(5, 8, 4, 2))  # Increase left margin for gene names
barplot(top10, 
        main = "Top 10 genes with most age-associated DMRs",
        xlab = "Number of DMRs",
        horiz = TRUE, 
        las = 1, 
        col = "steelblue",
        cex.names = 0.9)

dev.off()

# ====================================================
# 6. Figure 5: Summary table as image
# ====================================================
cat("Generating Figure 5: Summary table...\n")

# Create a summary table
summary_data <- data.frame(
  Metric = c("Total DMRs", "DMRs in genes", "PC1 variance", "PC2 age association (R²)",
             "Top gene 1", "Top gene 2", "Top gene 3"),
  Value = c("5,852", "80.4%", "18.5%", "21.2%",
            names(top10)[1], names(top10)[2], names(top10)[3])
)

# Save as CSV (useful for reference)
write.csv(summary_data, file.path(fig_dir, "05_summary_table.csv"), row.names = FALSE)

# Create a simple text file with key findings
sink(file.path(fig_dir, "05_key_findings.txt"))
cat("EPICv2 BRAIN TUMOR AGING ANALYSIS - KEY FINDINGS\n")
cat("==============================================\n\n")
cat("Total DMRs identified:", nrow(dmr_table), "\n")
cat("DMRs located in genes:", round(sum(dmr_table$overlapping_genes != "")/nrow(dmr_table)*100, 1), "%\n")
cat("\nTop 10 genes with most DMRs:\n")
print(top10)
cat("\nPC1 variance explained:", round(summary(pca)$importance[2,1]*100, 1), "%\n")
cat("PC2 age association (R²): 21.2%\n")
sink()

# ====================================================
# 7. Summary
# ====================================================
cat("\n========================================\n")
cat("FIGURES GENERATED SUCCESSFULLY\n")
cat("========================================\n")
cat("Files saved in:", fig_dir, "\n\n")
print(list.files(fig_dir))
# corrections after checking
# ====================================================
# Figure 1: Density plot (UPDATED VERSION)
# ====================================================
cat("Generating Figure 1: Density plot (updated)...\n")

png(file.path(fig_dir, "01_density_plot_updated.png"), width = 800, height = 600, res = 100)

plot(density(beta[,1], na.rm = TRUE), 
     col = 1, lwd = 2,
     main = "Beta values distribution across 16 samples\n(Each color represents one sample)",
     xlab = "Beta (methylation level)", 
     ylab = "Density",
     ylim = c(0, 5))

for(i in 2:ncol(beta)) {
  lines(density(beta[,i], na.rm = TRUE), col = i, lwd = 1.5)
}

# Nessuna legenda!

dev.off()

# ====================================================
# Figure 6: PC2 vs Age with regression line
# ====================================================
cat("Generating Figure 6: PC2 vs Age...\n")

png(file.path(fig_dir, "06_pc2_vs_age.png"), width = 800, height = 600, res = 100)

# Calcola R²
lm_pc2 <- lm(pca$x[,2] ~ sample_info$Age_numeric)
r2 <- round(summary(lm_pc2)$r.squared * 100, 1)

# Scatter plot
plot(sample_info$Age_numeric, pca$x[,2],
     xlab = "Age (years)", 
     ylab = "PC2",
     main = paste0("PC2 vs Age (R² = ", r2, "%)"),
     pch = 19, cex = 1.5,
     col = "steelblue")

# Aggiungi linea di regressione
abline(lm_pc2, col = "red", lwd = 2)

# Aggiungi etichette campioni (opzionale)
text(sample_info$Age_numeric, pca$x[,2], 
     labels = sample_info$Sample_Name, 
     pos = 3, cex = 0.7)

dev.off()

# ====================================================
# Figure 6: PC2 vs Age with regression line
# ====================================================
cat("Generating Figure 6: PC2 vs Age...\n")

png(file.path(fig_dir, "06_pc2_vs_age.png"), width = 800, height = 600, res = 100)

# Calcola R²
lm_pc2 <- lm(pca$x[,2] ~ sample_info$Age_numeric)
r2 <- round(summary(lm_pc2)$r.squared * 100, 1)

# Scatter plot
plot(sample_info$Age_numeric, pca$x[,2],
     xlab = "Age (years)", 
     ylab = "PC2",
     main = paste0("PC2 vs Age (R² = ", r2, "%)"),
     pch = 19, cex = 1.5,
     col = "steelblue",
     xlim = range(sample_info$Age_numeric) * c(0.95, 1.05))  # Leggero margine

# Aggiungi linea di regressione
abline(lm_pc2, col = "red", lwd = 2)

# Aggiungi numeri 1-16 invece dei codici GSM
text(sample_info$Age_numeric, pca$x[,2], 
     labels = 1:16, 
     pos = 3, cex = 0.9, col = "darkblue")

# Opzionale: aggiungi legenda con R²
legend("topleft", 
       legend = paste0("R² = ", r2, "%"), 
       bty = "n", cex = 1.2)

dev.off()
R.version.string
# ====================================================
# PCA colored by Tumor Type (Disease)
# ====================================================

# Load data
load("C:/Users/elega/OneDrive/Desktop/EPIC Project/data_processed/EPICv2_beta_values.RData")
sample_info <- read.csv("C:/Users/elega/OneDrive/Desktop/EPIC Project/metadata/Metadata_sample_sheet_epicv2_GPL33022.csv")

# Run PCA (if not already in memory)
pca <- prcomp(t(na.omit(beta)), center = TRUE, scale. = FALSE)

# Create a simplified tumor type grouping (if too many unique types)
# First, see what we have
table(sample_info$Disease)


# Installa dplyr (solo la prima volta)
install.packages("dplyr")

# Carica dplyr
library(dplyr)

# Poi riprova il codice
sample_info$Tumor_Group <- case_when(
  grepl("astrocytoma|glioma|glioblastoma", sample_info$Disease, ignore.case = TRUE) ~ "Glioma",
  grepl("medulloblastoma", sample_info$Disease, ignore.case = TRUE) ~ "Medulloblastoma",
  grepl("meningioma", sample_info$Disease, ignore.case = TRUE) ~ "Meningioma",
  TRUE ~ "Other"
)

# Verifica
table(sample_info$Tumor_Group)

# Crea dataframe per PCA
pca_df <- data.frame(
  PC1 = pca$x[,1],
  PC2 = pca$x[,2],
  Tumor_Group = sample_info$Tumor_Group,
  Sample = colnames(beta),
  Age = sample_info$Age_numeric
)

# Plot con ggplot2 (se installato)
if(require(ggplot2)) {
  ggplot(pca_df, aes(x = PC1, y = PC2, color = Tumor_Group)) +
    geom_point(size = 4) +
    geom_text(aes(label = 1:16), hjust = -0.3, vjust = 0.5, size = 3) +
    labs(
      title = "PCA colored by Tumor Group",
      x = paste0("PC1 (", round(summary(pca)$importance[2,1]*100, 1), "%)"),
      y = paste0("PC2 (", round(summary(pca)$importance[2,2]*100, 1), "%)")
    ) +
    theme_minimal()
  
  # Salva
  ggsave("C:/Users/elega/OneDrive/Desktop/EPIC Project/GitHub_Repo/figures/08_pca_tumor_grouped.png", 
         width = 10, height = 6, dpi = 300)
} else {
  # Versione base R
  colors <- as.numeric(factor(pca_df$Tumor_Group))
  plot(pca_df$PC1, pca_df$PC2, 
       col = colors, pch = 19, cex = 1.5,
       xlab = paste0("PC1 (", round(summary(pca)$importance[2,1]*100, 1), "%)"),
       ylab = paste0("PC2 (", round(summary(pca)$importance[2,2]*100, 1), "%)"),
       main = "PCA colored by Tumor Group")
  legend("topright", 
         legend = levels(factor(pca_df$Tumor_Group)), 
         col = 1:length(levels(factor(pca_df$Tumor_Group))), 
         pch = 19, cex = 1)
}
# ====================================================
# Calcolare R² per Tumor Type su PC1 e PC2
# ====================================================

# ANOVA per PC1 ~ Tumor_Group
aov_pc1 <- aov(pca$x[,1] ~ sample_info$Tumor_Group)
summary_pc1 <- summary(aov_pc1)

# ANOVA per PC2 ~ Tumor_Group  
aov_pc2 <- aov(pca$x[,2] ~ sample_info$Tumor_Group)
summary_pc2 <- summary(aov_pc2)

# Estrai R² (proporzione varianza spiegata)
r2_pc1 <- summary_pc1[[1]]$`Sum Sq`[1] / sum(summary_pc1[[1]]$`Sum Sq`)
r2_pc2 <- summary_pc2[[1]]$`Sum Sq`[1] / sum(summary_pc2[[1]]$`Sum Sq`)

# P-value
p_pc1 <- summary_pc1[[1]]$`Pr(>F)`[1]
p_pc2 <- summary_pc2[[1]]$`Pr(>F)`[1]

# Risultati
cat("=== TUMOR TYPE ASSOCIATION ===\n")
cat("PC1 ~ Tumor Type:\n")
cat("  R² =", round(r2_pc1 * 100, 1), "%\n")
cat("  p-value =", format(p_pc1, scientific = TRUE, digits = 3), "\n\n")

cat("PC2 ~ Tumor Type:\n")
cat("  R² =", round(r2_pc2 * 100, 1), "%\n")
cat("  p-value =", format(p_pc2, scientific = TRUE, digits = 3), "\n")

# Confronto con Age
cat("\n=== CONFRONTO CON AGE ===\n")
cat("PC2 ~ Age: R² = 21.2%, p < 0.05\n")


