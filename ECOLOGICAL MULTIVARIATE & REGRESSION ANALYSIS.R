# ==============================================================================
# ECOLOGICAL MULTIVARIATE & REGRESSION ANALYSIS (CERRADO)
# Methods: NMDS, Environmental PCA, Simple & Partial Mantel Tests, Poisson GLMs
# ==============================================================================

# ------------------------------------------------------------------------------
# 1. PACKAGE INSTALLATION & LOADING
# ------------------------------------------------------------------------------
required_packages <- c("vegan", "geosphere", "tidyverse", "ggrepel", "gridExtra")
missing_packages  <- required_packages[!(required_packages %in% installed.packages()[, "Package"])]

if (length(missing_packages) > 0) {
  install.packages(missing_packages)
}

library(vegan)
library(geosphere)
library(tidyverse)
library(ggrepel)
library(gridExtra)

# Set seed for reproducibility in permutation tests and NMDS
set.seed(42)

# ------------------------------------------------------------------------------
# 2. DATA IMPORT & ALIGNMENT
# ------------------------------------------------------------------------------
cat("\n=== 1. LOADING AND ALIGNING DATA ===\n")

floristic     <- read.csv("areas cerrado - Planilha1.csv", row.names = 1, check.names = FALSE)
environmental <- read.csv("areas_cerrado_com_variaveis_GDM.csv", row.names = 1)

# Keep only overlapping sites and align row order
common_sites  <- intersect(rownames(floristic), rownames(environmental))
floristic     <- floristic[common_sites, ]
environmental <- environmental[common_sites, ]

# Separate geographic coordinates from environmental predictors
coordinates <- environmental[, c("longitude", "latitude")]
env_vars    <- environmental[, !names(environmental) %in% c("latitude", "longitude")]

cat(paste("[OK] Data aligned:", nrow(floristic), "sites and", ncol(floristic), "species.\n"))

# ------------------------------------------------------------------------------
# 3. FLORISTIC NMDS ORDINATION
# ------------------------------------------------------------------------------
cat("\n=== 2. RUNNING FLORISTIC NMDS ===\n")

# Calculate Jaccard dissimilarity matrix (binary presence/absence)
dist_floristic <- vegdist(floristic, method = "jaccard", binary = TRUE)

# Run 2D NMDS
nmds_res    <- metaMDS(dist_floristic, k = 2, trymax = 100, trace = FALSE)
scores_nmds <- as.data.frame(scores(nmds_res, display = "sites"))
scores_nmds$Area <- rownames(scores_nmds)

# Save NMDS plot
png("nmds_floristic_plot.png", width = 4000, height = 3000, res = 500)

plot_nmds <- ggplot(scores_nmds, aes(x = NMDS1, y = NMDS2)) +
  geom_point(color = "#2E8B57", size = 4, alpha = 0.8) +
  geom_text_repel(
    aes(label = Area), size = 3.5, fontface = "bold", color = "grey40",
    box.padding = 0.5, point.padding = 0.3, max.overlaps = Inf
  ) +
  annotate(
    "text", x = -Inf, y = Inf,
    label = paste("Stress =", round(nmds_res$stress, 4)),
    hjust = -0.1, vjust = 1.5, size = 4.5, fontface = "italic", color = "black"
  ) +
  labs(title = "", x = "NMDS1", y = "NMDS2") +
  theme_bw(base_size = 14) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5),
    axis.text  = element_text(color = "black")
  )

print(plot_nmds)
dev.off()
cat("[SUCCESS] NMDS plot saved as 'nmds_floristic_plot.png'.\n")

# ------------------------------------------------------------------------------
# 4. ENVIRONMENTAL PCA
# ------------------------------------------------------------------------------
cat("\n=== 3. RUNNING ENVIRONMENTAL PCA ===\n")

env_scaled <- scale(env_vars)
pca_res    <- prcomp(env_scaled, scale. = TRUE)
scores_pca <- as.data.frame(pca_res$x[, 1:2])
colnames(scores_pca) <- c("PC1", "PC2")
scores_pca$Area <- rownames(scores_pca)

var_explained <- (pca_res$sdev^2) / sum(pca_res$sdev^2) * 100

# Save PCA plot
png("pca_environmental_plot.png", width = 4000, height = 3000, res = 500)

plot_pca <- ggplot(scores_pca, aes(x = PC1, y = PC2)) +
  geom_point(color = "#4682B4", size = 4, alpha = 0.8) +
  geom_text_repel(
    aes(label = Area), size = 3.5, fontface = "bold",
    box.padding = 0.5, point.padding = 0.3, max.overlaps = Inf
  ) +
  labs(
    title = "PCA Ordination of Environmental Variables (Soil & Climate)",
    x = paste0("Principal Component 1 (", round(var_explained[1], 2), "%)"),
    y = paste0("Principal Component 2 (", round(var_explained[2], 2), "%)")
  ) +
  theme_bw(base_size = 14) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5),
    axis.text  = element_text(color = "black")
  )

print(plot_pca)
dev.off()
cat("[SUCCESS] PCA plot saved as 'pca_environmental_plot.png'.\n")

# ------------------------------------------------------------------------------
# 5. MATRIX CORRELATIONS: SIMPLE & PARTIAL MANTEL TESTS
# ------------------------------------------------------------------------------
cat("\n=== 4. RUNNING SIMPLE & PARTIAL MANTEL TESTS ===\n")

# Compute Euclidean environmental distance and Haversine geographic distance
dist_environmental <- dist(env_scaled, method = "euclidean")
geo_matrix         <- distm(coordinates, fun = distHaversine)
dist_geographic    <- as.dist(geo_matrix)

# A. Simple Mantel: Floristic vs. Environmental
mantel_env <- mantel(dist_floristic, dist_environmental, method = "pearson", permutations = 999)

# B. Simple Mantel: Floristic vs. Geographic
mantel_geo <- mantel(dist_floristic, dist_geographic, method = "pearson", permutations = 999)

# C. Partial Mantel: Floristic vs. Environmental (controlling for Space)
mantel_partial_env <- mantel.partial(dist_floristic, dist_environmental, dist_geographic, method = "pearson", permutations = 999)

# D. Partial Mantel: Floristic vs. Geographic (controlling for Environment)
mantel_partial_geo <- mantel.partial(dist_floristic, dist_geographic, dist_environmental, method = "pearson", permutations = 999)

# Prepare distance vectors for plotting (converting meters to km)
df_distances <- data.frame(
  Floristic_Dist     = as.vector(dist_floristic),
  Environmental_Dist = as.vector(dist_environmental),
  Geographic_Dist    = as.vector(dist_geographic) / 1000
)

# Dynamic labels for plot annotations
label_env <- paste0("r = ", round(mantel_env$statistic, 4), "\np = ", mantel_env$signif)
label_geo <- paste0("r = ", round(mantel_geo$statistic, 4), "\np = ", mantel_geo$signif)

# Plot A: Floristic vs. Environmental Distance
plot_mantel_env <- ggplot(df_distances, aes(x = Environmental_Dist, y = Floristic_Dist)) +
  geom_point(color = "#2E8B57", alpha = 0.6, size = 2.5) +
  geom_smooth(method = "lm", color = "#114B28", se = TRUE, linewidth = 1.2) +
  annotate(
    "text", x = Inf, y = -Inf, label = label_env,
    hjust = 1.1, vjust = -0.5, size = 4.5, fontface = "italic", color = "black"
  ) +
  labs(title = "", x = "Environmental Distance", y = "Floristic Distance") +
  theme_bw(base_size = 13) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5),
    axis.text  = element_text(color = "black")
  )

# Plot B: Floristic vs. Geographic Distance
plot_mantel_geo <- ggplot(df_distances, aes(x = Geographic_Dist, y = Floristic_Dist)) +
  geom_point(color = "#2E8B57", alpha = 0.6, size = 2.5) +
  geom_smooth(method = "lm", color = "#114B28", se = TRUE, linewidth = 1.2) +
  annotate(
    "text", x = Inf, y = -Inf, label = label_geo,
    hjust = 1.1, vjust = -0.5, size = 4.5, fontface = "italic", color = "black"
  ) +
  labs(title = "", x = "Geographic Distance (km)", y = "Floristic Distance") +
  theme_bw(base_size = 13) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5),
    axis.text  = element_text(color = "black")
  )

# Save combined Mantel plot
png("mantel_plots_combined.png", width = 3000, height = 1500, res = 300)
grid.arrange(plot_mantel_env, plot_mantel_geo, ncol = 2)
dev.off()
cat("[SUCCESS] Combined Mantel plot saved as 'mantel_plots_combined.png'.\n")

# ------------------------------------------------------------------------------
# 6. GENERALIZED LINEAR MODELS (POISSON GLM) - RICHNESS vs. ENVIRONMENT
# ------------------------------------------------------------------------------
cat("\n=== 5. RUNNING POISSON GLMs FOR SPECIES RICHNESS ===\n")

richness <- rowSums(floristic > 0)
glm_data <- data.frame(Richness = richness, env_vars)

# Reshape data to long format for faceted plotting
glm_data_long <- glm_data %>%
  pivot_longer(
    cols = -Richness,
    names_to = "Variable",
    values_to = "Value"
  )

# Clean English labels for environmental variables
clean_labels_en <- c(
  "Coefficient_of_variation_EVI_01_05_1km_uint16" = "EVI variation",
  "et0_v31_yr"                                    = "Potential evapotranspiration",
  "fogo"                                          = "Fire frequency",
  "soilgrids.isric_myRegion_clay_5.15cm_mean"     = "Clay content",
  "soilgrids.isric_myRegion_phh2o_5.15cm_mean"    = "Soil water pH",
  "soilgrids.isric_myRegion_sand_5.15cm_mean"    = "Sand content",
  "soilgrids.isric_myRegion_silt_5.15cm_mean"    = "Silt content",
  "terrain_ruggedness_index_1KMmd_GMTEDmd"        = "Terrain ruggedness index",
  "wc2.1_30s_bio_1"                               = "Bio 1",
  "wc2.1_30s_bio_14"                              = "Bio 14",
  "wc2.1_30s_bio_18"                              = "Bio 18",
  "wc2.1_30s_bio_3"                               = "Bio 3",
  "wc2.1_30s_bio_4"                               = "Bio 4",
  "wc2.1_30s_bio_8"                               = "Bio 8",
  "wc2.1_30s_bio_9"                               = "Bio 9"
)

glm_data_long$Variable_Label <- clean_labels_en[glm_data_long$Variable]

# Save multi-panel GLM curves plot
png("all_glms_poisson_curves.png", width = 4000, height = 3200, res = 500)

plot_glms <- ggplot(glm_data_long, aes(x = Value, y = Richness)) +
  geom_point(color = "#2E8B57", alpha = 0.6, size = 2) +
  geom_smooth(
    method = "glm", method.args = list(family = "poisson"),
    color = "#114B28", fill = "#2E8B57", alpha = 0.15, linewidth = 1
  ) +
  facet_wrap(~ Variable_Label, ncol = 5, scales = "free_x") +
  labs(
    title = "",
    subtitle = "",
    x = "Environmental Predictor Values",
    y = "Species Richness"
  ) +
  theme_bw(base_size = 10) +
  theme(
    strip.background = element_rect(fill = "#f0f0f0", color = "black"),
    strip.text       = element_text(face = "bold", size = 7, color = "black"),
    axis.text        = element_text(color = "black"),
    panel.grid.minor = element_blank()
  )

print(plot_glms)
dev.off()
cat("[SUCCESS] Multi-panel GLM curves saved as 'all_glms_poisson_curves.png'.\n")

# ------------------------------------------------------------------------------
# 7. EXPORT STATISTICAL SUMMARY REPORT (.TXT)
# ------------------------------------------------------------------------------
cat("\n=== 6. GENERATING STATISTICAL REPORT ===\n")

sink("ecological_analysis_report.txt")

cat("=====================================================================\n")
cat("             ECOLOGICAL ANALYSIS SUMMARY REPORT (CERRADO)            \n")
cat("          (NMDS, Environmental PCA, Mantel Tests & Poisson GLMs)     \n")
cat("=====================================================================\n\n")

cat(paste("Total Sites Analyzed:", nrow(floristic), "\n"))
cat(paste("Total Species Recorded:", ncol(floristic), "\n\n"))

# 1. NMDS
cat("---------------------------------------------------------------------\n")
cat("1. FLORISTIC NMDS ORDINATION\n")
cat("---------------------------------------------------------------------\n")
cat("Dissimilarity Metric: Jaccard Index (binary presence/absence)\n")
cat(paste("NMDS Stress (2 Dimensions):", round(nmds_res$stress, 4), "\n\n"))

# 2. PCA
cat("---------------------------------------------------------------------\n")
cat("2. ENVIRONMENTAL PCA EXPLAINED VARIANCE\n")
cat("---------------------------------------------------------------------\n")
cat(paste0("Variance explained by Axis 1 (PC1): ", round(var_explained[1], 2), "%\n"))
cat(paste0("Variance explained by Axis 2 (PC2): ", round(var_explained[2], 2), "%\n\n"))

# 3. Mantel Tests
cat("---------------------------------------------------------------------\n")
cat("3. MATRIX CORRELATION ANALYSIS (MANTEL TESTS)\n")
cat("---------------------------------------------------------------------\n")
cat("A. Simple Mantel Tests:\n")
cat(paste0("  - Floristic vs. Environment: r = ", round(mantel_env$statistic, 4), " | p-value = ", mantel_env$signif, "\n"))
cat(paste0("  - Floristic vs. Geography:   r = ", round(mantel_geo$statistic, 4), " | p-value = ", mantel_geo$signif, "\n\n"))

cat("B. Partial Mantel Tests (Isolated Effects):\n")
cat(paste0("  - Floristic vs. Environment (Controlling Space):       r = ", round(mantel_partial_env$statistic, 4), " | p-value = ", mantel_partial_env$signif, "\n"))
cat(paste0("  - Floristic vs. Geography   (Controlling Environment): r = ", round(mantel_partial_geo$statistic, 4), " | p-value = ", mantel_partial_geo$signif, "\n\n"))

cat("Ecological Interpretation Guidelines:\n")
cat(" - Significant r (Floristic vs. Environment | Space): Evidence of Niche Selection / Environmental Filtering.\n")
cat(" - Significant r (Floristic vs. Geography | Environment): Evidence of Dispersal Limitation.\n\n")

# 4. Poisson GLMs
cat("---------------------------------------------------------------------\n")
cat("4. GENERALIZED LINEAR MODELS (POISSON GLM) - RICHNESS vs. ENVIRONMENT\n")
cat("---------------------------------------------------------------------\n")

glm_summary_df <- data.frame(
  Variable    = character(),
  Clean_Name  = character(),
  Coefficient = numeric(),
  P_Value     = numeric(),
  Significant = character()
)

for (var in names(env_vars)) {
  current_formula <- as.formula(paste("Richness ~", var))
  model_fit       <- glm(current_formula, data = glm_data, family = poisson)
  model_summary   <- summary(model_fit)
  
  coef_val <- round(model_summary$coefficients[2, 1], 5)
  p_val    <- round(model_summary$coefficients[2, 4], 5)
  sig_flag <- ifelse(p_val < 0.05, "Yes *", "No")
  
  glm_summary_df <- rbind(glm_summary_df, data.frame(
    Variable    = var,
    Clean_Name  = clean_labels_en[var],
    Coefficient = coef_val,
    P_Value     = p_val,
    Significant = sig_flag
  ))
}

print(glm_summary_df, row.names = FALSE)

cat("\n=====================================================================\n")
cat("END OF REPORT\n")
cat("=====================================================================\n")

sink()

cat("[SUCCESS] Report generated as 'ecological_analysis_report.txt'!\n")
