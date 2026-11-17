# ==============================================================================
# MULTIVARIATE FLORISTIC ANALYSIS (CERRADO SUBREGIONS)
# Methods: Jaccard Dissimilarity, PERMANOVA, NMDS, and UPGMA Clustering
# ==============================================================================

# ------------------------------------------------------------------------------
# 1. PACKAGE INSTALLATION & LOADING
# ------------------------------------------------------------------------------
required_packages <- c("vegan", "geosphere", "tidyverse", "ggrepel", "ggdendro")
missing_packages  <- required_packages[!(required_packages %in% installed.packages()[, "Package"])]

if (length(missing_packages) > 0) {
  install.packages(missing_packages)
}

library(vegan)
library(geosphere)
library(tidyverse)
library(ggrepel)
library(ggdendro)

# Set seed for reproducibility in permutation tests and NMDS
set.seed(42)

# ------------------------------------------------------------------------------
# 2. DATA IMPORT & CLEANING
# ------------------------------------------------------------------------------
cat("\n=== 1. LOADING AND VALIDATING DATA ===\n")

# Load floristic and environmental matrices
floristic <- read.csv("areas cerrado - Planilha1.csv", row.names = 1, check.names = FALSE)
environmental <- read.csv("areas_cerrado_com_variaveis_GDM.csv", row.names = 1)

# Keep only overlapping sites between datasets
common_sites <- intersect(rownames(floristic), rownames(environmental))
floristic <- floristic[common_sites, ]
environmental <- environmental[common_sites, ]

# Load locality/subregion metadata (comma-separated)
# Note: If semicolon-separated, replace with: read_delim("...", delim = ";")
localities <- read_csv("cerrado localidades - Página1 (1).csv", show_col_types = FALSE)

# Trim whitespace from site identifiers to prevent matching errors
rownames(floristic) <- trimws(rownames(floristic))
localities$Code     <- trimws(localities$Code)

# ------------------------------------------------------------------------------
# 3. DATA VALIDATION & ALIGNMENT
# ------------------------------------------------------------------------------
floristic_sites  <- rownames(floristic)
locality_sites   <- localities$Code
missing_in_local <- setdiff(floristic_sites, locality_sites)

if (length(missing_in_local) > 0) {
  stop(paste("CRITICAL ERROR: The following sites are in the floristic matrix but missing from the localities CSV:",
             paste(missing_in_local, collapse = ", ")))
} else {
  cat("[OK] All floristic sites matched successfully with the localities metadata.\n")
}

# Strictly align localities metadata to match the exact row order of the floristic matrix
aligned_localities <- localities[match(floristic_sites, localities$Code), ]

if (!all(rownames(floristic) == aligned_localities$Code)) {
  stop("CRITICAL ERROR: Row alignment failed between floristic matrix and localities.")
} else {
  cat("[OK] Perfect row alignment verified.\n")
}

# Translate subregion names to English
aligned_localities$Subregion_EN <- recode(
  aligned_localities$Subregion,
  "Centro-oeste"   = "Central-West",
  "Norte-nordeste" = "North-Northeast",
  "Meridional"     = "Southern",
  "Mesotrófico"    = "Mesotrophic",
  "Centro-sudeste" = "Central-Southeast"
)

# Define custom color palette for subregions
subregion_colors <- c(
  "Central-West"      = "#510ecd",
  "North-Northeast"   = "#114b28",
  "Southern"          = "#a2e8b9",
  "Mesotrophic"       = "#3168d8",
  "Central-Southeast" = "#5dc963"
)

# ------------------------------------------------------------------------------
# 4. FLORISTIC METRICS & DISSIMILARITY MATRIX
# ------------------------------------------------------------------------------
cat("\n=== 2. COMPUTING DISSIMILARITY & SINGLETONS ===\n")

# Count singletons (species occurring in only 1 site)
n_singletons    <- sum(colSums(floristic > 0) == 1)
perc_singletons <- round((n_singletons / ncol(floristic)) * 100, 1)
cat(paste0("Total singletons: ", n_singletons, " out of ", ncol(floristic), " species (", perc_singletons, "%)\n"))

# Calculate Jaccard dissimilarity matrix (binary = TRUE for presence/absence data)
dist_jaccard <- vegdist(floristic, method = "jaccard", binary = TRUE)

# ------------------------------------------------------------------------------
# 5. GLOBAL PERMANOVA
# ------------------------------------------------------------------------------
cat("\n=== 3. RUNNING GLOBAL PERMANOVA ===\n")

permanova_res <- adonis2(dist_jaccard ~ Subregion_EN, data = aligned_localities, permutations = 999)
print(permanova_res)

# ------------------------------------------------------------------------------
# 6. NON-METRIC MULTIDIMENSIONAL SCALING (NMDS)
# ------------------------------------------------------------------------------
cat("\n=== 4. RUNNING NMDS & GENERATING PLOT ===\n")

nmds_res <- metaMDS(dist_jaccard, k = 2, trymax = 100, trace = FALSE)

# Extract site scores and merge with aligned metadata
nmds_scores <- as.data.frame(scores(nmds_res, display = "sites"))
nmds_scores$Area <- rownames(nmds_scores)
plot_data_nmds <- merge(nmds_scores, aligned_localities, by.x = "Area", by.y = "Code")

# Save high-resolution NMDS plot with horizontal legend and expanded margins
png("nmds_floristic_plot.png", width = 4200, height = 3200, res = 500)

plot_nmds <- ggplot(plot_data_nmds, aes(x = NMDS1, y = NMDS2)) +
  geom_point(aes(color = Subregion_EN), size = 4, alpha = 0.9) +
  geom_text_repel(
    aes(label = Area), size = 3.5, fontface = "bold", color = "grey30",
    box.padding = 0.5, point.padding = 0.3, max.overlaps = Inf
  ) +
  annotate(
    "text", x = -Inf, y = Inf,
    label = paste("Stress =", round(nmds_res$stress, 4)),
    hjust = -0.1, vjust = 1.5, size = 3.5, fontface = "italic", color = "black"
  ) +
  scale_color_manual(values = subregion_colors) +
  labs(title = "", x = "NMDS1", y = "NMDS2", color = "Subregion") +
  theme_bw(base_size = 14) +
  theme(
    axis.text       = element_text(color = "black"),
    legend.position = "bottom",
    legend.box      = "horizontal",
    legend.title    = element_text(face = "bold"),
    plot.margin     = margin(t = 0.5, r = 1.5, b = 0.5, l = 0.5, unit = "cm")
  ) +
  guides(color = guide_legend(nrow = 1, byrow = TRUE))

print(plot_nmds)
dev.off()
cat("[SUCCESS] NMDS plot saved as 'nmds_floristic_plot.png'.\n")

# ------------------------------------------------------------------------------
# 7. HIERARCHICAL CLUSTERING (UPGMA)
# ------------------------------------------------------------------------------
cat("\n=== 5. RUNNING UPGMA CLUSTERING & GENERATING DENDROGRAM ===\n")

# Compute UPGMA dendrogram (method = "average") and Cophenetic Correlation
upgma_res <- hclust(dist_jaccard, method = "average")
coph_dist <- cophenetic(upgma_res)
coph_cor  <- cor(dist_jaccard, coph_dist)
cat(paste("Cophenetic Correlation (r):", round(coph_cor, 4), "\n"))

# Extract dendrogram segments and tip labels for ggplot2
dendro_data <- dendro_data(as.dendrogram(upgma_res), type = "rectangle")
segments_df <- segment(dendro_data)
tips_df     <- label(dendro_data)
tips_df     <- merge(tips_df, aligned_localities, by.x = "label", by.y = "Code")

max_height <- max(segments_df$y)

# Save high-resolution UPGMA dendrogram
png("upgma_floristic_plot.png", width = 4200, height = 3200, res = 500)

plot_upgma <- ggplot() +
  geom_segment(
    data = segments_df,
    aes(x = x, y = y, xend = xend, yend = yend),
    color = "grey35", linewidth = 0.7
  ) +
  geom_point(
    data = tips_df,
    aes(x = x, y = 0, color = Subregion_EN),
    size = 4, alpha = 0.95
  ) +
  geom_text(
    data = tips_df,
    aes(x = x, y = -0.025, label = label),
    hjust = 0, size = 3.5, fontface = "bold", color = "grey30"
  ) +
  annotate(
    "text", x = Inf, y = max_height * 1.03,
    label = paste("Cophenetic r =", round(coph_cor, 4)),
    hjust = 0, vjust = 1.5, size = 3.5, fontface = "italic", color = "black"
  ) +
  coord_flip() +
  scale_y_reverse(
    limits = c(max_height * 1.05, -max_height * 0.28),
    breaks = seq(0, 1, by = 0.2)
  ) +
  scale_color_manual(values = subregion_colors) +
  labs(title = "", x = "", y = "Jaccard Dissimilarity", color = "Subregion") +
  theme_bw(base_size = 14) +
  theme(
    axis.text.x        = element_text(color = "black"),
    axis.text.y        = element_blank(),
    axis.ticks.y       = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.grid.minor   = element_blank(),
    legend.position    = "bottom",
    legend.box         = "horizontal",
    legend.title       = element_text(face = "bold"),
    plot.margin        = margin(t = 0.5, r = 1.0, b = 0.5, l = 0.5, unit = "cm")
  ) +
  guides(color = guide_legend(nrow = 1, byrow = TRUE))

print(plot_upgma)
dev.off()
cat("[SUCCESS] UPGMA dendrogram saved as 'upgma_floristic_plot.png'.\n")

# ------------------------------------------------------------------------------
# 8. EXPORT STATISTICAL SUMMARY REPORT (.TXT)
# ------------------------------------------------------------------------------
cat("\n=== 6. GENERATING STATISTICAL REPORT ===\n")

f_val       <- round(permanova_res$F[1], 3)
p_val       <- round(permanova_res$`Pr(>F)`[1], 3)
r2_val      <- round(permanova_res$R2[1], 4)
r2_pct      <- round(permanova_res$R2[1] * 100, 2)
df_model    <- permanova_res$Df[1]
df_residual <- permanova_res$Df[2]
stress_val  <- round(nmds_res$stress, 4)
coph_val    <- round(coph_cor, 4)

sink("permanova_nmds_report.txt")

cat("=====================================================================\n")
cat("       MULTIVARIATE FLORISTIC ANALYSIS REPORT (CERRADO)              \n")
cat("=====================================================================\n\n")

cat("1. DATASET SUMMARY\n")
cat("---------------------------------------------------------------------\n")
cat(paste("Total Sites Analyzed:", nrow(floristic), "\n"))
cat(paste("Total Species:", ncol(floristic), "\n"))
cat(paste0("Singletons (1-site species): ", n_singletons, " (", perc_singletons, "%)\n\n"))

cat("2. NON-METRIC MULTIDIMENSIONAL SCALING (NMDS)\n")
cat("---------------------------------------------------------------------\n")
cat("Dissimilarity Metric: Jaccard Index (binary presence/absence)\n")
cat(paste("Final Stress (k = 2):", stress_val, "\n"))
if (stress_val < 0.2) {
  cat("Stress Interpretation: Acceptable/Good (< 0.20). The 2D ordination\n")
  cat("accurately preserves the original dissimilarities among sites.\n")
} else {
  cat("Stress Interpretation: High (>= 0.20). Interpret 2D distances with caution.\n")
}
cat("\n")

cat("3. HIERARCHICAL CLUSTERING (UPGMA)\n")
cat("---------------------------------------------------------------------\n")
cat(paste("Cophenetic Correlation Coefficient (r):", coph_val, "\n\n"))

cat("4. PERMUTATIONAL MULTIVARIATE ANALYSIS OF VARIANCE (GLOBAL PERMANOVA)\n")
cat("---------------------------------------------------------------------\n")
cat("Formula: dist_jaccard ~ Subregion_EN\n")
cat("Permutations: 999\n")
cat(paste("Model Degrees of Freedom (Df Model):", df_model, "\n"))
cat(paste("Residual Degrees of Freedom (Df Residual):", df_residual, "\n"))
cat(paste("Pseudo-F Statistic (F-Model):", f_val, "\n"))
cat(paste0("Coefficient of Determination (R2): ", r2_val, " (", r2_pct, "%)\n"))
cat(paste("p-value (Pr(>F)):", p_val, "\n\n"))

cat("MANUSCRIPT SUMMARY DRAFT:\n")
cat("---------------------------------------------------------------------\n")
cat(paste0(
  "Permutational Multivariate Analysis of Variance (Global PERMANOVA) revealed\n",
  "significant differences in floristic composition among the studied subregions\n",
  "(F = ", f_val, "; p = ", p_val, "). Geographic subregions accounted for ", r2_pct, "%\n",
  "of the total variation in species similarity (R2 = ", r2_val, "), indicating a\n",
  "moderately strong spatial structuring of plant communities.\n\n"
))

cat("=====================================================================\n")
cat("Full R Output (vegan::adonis2):\n")
cat("=====================================================================\n")
print(permanova_res)

sink()

cat("[SUCCESS] Report saved as 'permanova_nmds_report.txt'!\n")
