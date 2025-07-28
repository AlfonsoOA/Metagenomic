# ====================
# 0. Libraries and Setup
# ====================
packages <- c("tidyverse", "vegan", "ggpubr", "reshape2", "pheatmap",
              "VennDiagram", "ape", "igraph", "psych", "cluster", "dendextend", "Hmisc") # Añadido Hmisc explícitamente

invisible(lapply(packages, function(x) if (!require(x, character.only = TRUE)) install.packages(x)))
invisible(lapply(packages, library, character.only = TRUE))

# ===========================
# 1. Load and Prepare Data
# ===========================
base_path <- "C:/Users/Alfonso/Desktop/AlfonsoOA_MSI/1Investigacion/5BiodegradacionPlasticos/11MetagenomaTenebrio/0_Metagenome_Figure2"
file <- file.path(base_path, "1_EC_relative_abundance_table_tenebrio_molitor.csv")
out_dir <- file.path(base_path, "Metagenome_Results") # Renamed output directory
dir.create(out_dir, showWarnings = FALSE)

# Load data without factors to avoid issues
data_raw <- read.csv(file, sep = "\t", stringsAsFactors = FALSE)

# Set first column as row names
rownames(data_raw) <- data_raw$gene

# Select abundance columns (columns 5 to 16 = 3 conditions * 4 replicates)
abund_cols <- colnames(data_raw)[5:16]

# Extract abundance matrix: rows = genes, columns = samples
data_abund <- data_raw[, abund_cols]

# Convert to numeric matrix for analysis
data_abund <- as.matrix(data_abund)
mode(data_abund) <- "numeric"

# Define condition factor by extracting part before underscore
cond <- factor(gsub("_\\dK", "", abund_cols), levels = c("C", "LD", "PET"))

# Configuration for saving high-quality images
dpi <- 600
width_in <- 7
height_in <- 5

# ==================
# 2. Alpha Diversity
# ==================
alpha_div <- data.frame(
  sample = colnames(data_abund),
  condition = cond,
  Shannon = vegan::diversity(t(data_abund), index = "shannon"),
  Simpson = vegan::diversity(t(data_abund), index = "simpson")
)

write.csv(alpha_div, file.path(out_dir, "alpha_diversity.csv"), row.names = FALSE)

# Function to create and save alpha diversity plot
plot_alpha_diversity <- function(data, index_name, out_dir) {
  p <- ggplot(data, aes(x = condition, y = .data[[index_name]], fill = condition)) +
    geom_boxplot(outlier.size = 1) +
    geom_jitter(width = 0.2, alpha = 0.7, size = 1.5) +
    labs(
      title = paste(index_name, "Diversity by Condition"),
      x = "Condition",
      y = paste(index_name, "Diversity Index")
    ) +
    theme_minimal() +
    theme(legend.position = "none", text = element_text(size = 14))
  
  # Save in TIFF
  tiff_filename <- file.path(out_dir, paste0("alpha_diversity_", index_name, ".tiff"))
  tiff(tiff_filename, width = width_in, height = height_in, units = "in", res = dpi)
  print(p)
  dev.off()
  
  # Save in PNG
  png_filename <- file.path(out_dir, paste0("alpha_diversity_", index_name, ".png"))
  png(png_filename, width = width_in, height = height_in, units = "in", res = dpi)
  print(p)
  dev.off()
}

# Assuming 'alpha_div' and 'out_dir' are already defined
plot_alpha_diversity(alpha_div, "Shannon", out_dir)
plot_alpha_diversity(alpha_div, "Simpson", out_dir)

# =====================
# 3. Sample Rarefaction (simulada a partir de datos TMM)
# =====================

# Leer archivo sin fijar aún row.names
counts_raw <- read.table(file.path(base_path, "Abundance_Table.txt"), sep = "\t", header = TRUE, quote = "", check.names = FALSE)

# Verificar estructura
# str(counts_raw)  # Puedes descomentar para inspeccionar

# Usar columna 'gene' como identificador
rownames(counts_raw) <- counts_raw$gene

# Seleccionar solo columnas de abundancia
abund_cols <- grep("^C_\\dK$|^LD_\\dK$|^PET_\\dK$", colnames(counts_raw), value = TRUE)
raw_counts <- counts_raw[, abund_cols]

# Redondear (simulación de conteos enteros)
raw_counts <- round(raw_counts)

# Transponer para que filas = muestras, columnas = genes
data_for_rare <- t(raw_counts)

# Comprobar que no hay negativos
stopifnot(all(data_for_rare >= 0))

# Definir vector de condiciones
cond <- factor(gsub("_\\dK", "", rownames(data_for_rare)), levels = c("C", "LD", "PET"))

# Colores por condición
cols_cond <- c(C = "blue", LD = "red", PET = "green")

# Valor mínimo de lecturas totales por muestra
raremax <- min(rowSums(data_for_rare))

# Guardar curva de rarefacción en TIFF
tiff(filename = file.path(out_dir, "Rarefaction_curve.tiff"), width = width_in, height = height_in, units = "in", res = dpi)
rarecurve(data_for_rare, step = 100, sample = raremax, label = FALSE,
          col = cols_cond[cond], lwd = 1.5,
          xlab = "Reads per sample", ylab = "Genes detected", main = "Rarefaction Curves (Simulated from TMM)")
legend("bottomright", legend = levels(cond), col = cols_cond, lwd = 2, bty = "n")
dev.off()

# Guardar curva de rarefacción en PNG
png(filename = file.path(out_dir, "Rarefaction_curve.png"), width = width_in, height = height_in, units = "in", res = dpi)
rarecurve(data_for_rare, step = 100, sample = raremax, label = FALSE,
          col = cols_cond[cond], lwd = 1.5,
          xlab = "Reads per sample", ylab = "Genes detected", main = "Rarefaction Curves (Simulated from TMM)")
legend("bottomright", legend = levels(cond), col = cols_cond, lwd = 2, bty = "n")
dev.off()

# ===================
# 4. Beta Diversity (Ordination and Statistical Tests)
# ===================

# Calculate Bray-Curtis distance matrix between samples
dist_bc <- vegdist(t(data_abund), method = "bray")

# --- 4.1. Principal Coordinate Analysis (PCoA) and PERMANOVA
pcoa_res <- ape::pcoa(dist_bc)

# Prepare dataframe for plot
pcoa_df <- data.frame(
  Sample = colnames(data_abund),
  Condition = cond,
  Axis1 = pcoa_res$vectors[,1],
  Axis2 = pcoa_res$vectors[,2]
)

# Save PCoA table
write.csv(pcoa_df, file.path(out_dir, "PCoA_coordinates.csv"), row.names = FALSE)

# PERMANOVA for statistical test of differences (associated with PCoA)
adonis_res <- vegan::adonis2(dist_bc ~ cond, permutations = 999)
write.csv(as.data.frame(adonis_res), file.path(out_dir, "PERMANOVA_results.csv"))

# Extract PERMANOVA R2 and p-value for plot annotation
permanova_r2 <- round(adonis_res$R2[1], 3)
permanova_p <- format.pval(adonis_res$`Pr(>F)`[1], digits = 3) # Format p-value for display

# PCoA Plot with ggplot2, including PERMANOVA results and connected points
p_pcoa <- ggplot(pcoa_df, aes(x = Axis1, y = Axis2, color = Condition, group = Condition)) + # Added group = Condition
  geom_point(size = 4, aes(shape = Condition)) + # Added shape aesthetic
  stat_ellipse(level = 0.95) + # Add 95% confidence ellipses
  geom_line(aes(group = Condition), alpha = 0.3) + # Conectar puntos dentro de cada grupo
  theme_minimal() +
  labs(title = "PCoA of Gene Abundance Profiles",
       x = paste0("PCoA1 (", round(pcoa_res$values$Relative_eig[1]*100, 2), "%)"),
       y = paste0("PCoA2 (", round(pcoa_res$values$Relative_eig[2]*100, 2), "%)"),
       subtitle = paste0("PERMANOVA: R² = ", permanova_r2, ", p = ", permanova_p)) +
  theme(plot.subtitle = element_text(hjust = 0.5, size = 12)) # Center subtitle

# Save PCoA plot in PNG
png(filename = file.path(out_dir, "PCoA_plot.png"), width = width_in, height = height_in, units = "in", res = dpi)
print(p_pcoa)
dev.off()

# Save PCoA plot in TIFF
tiff(filename = file.path(out_dir, "PCoA_plot.tiff"), width = width_in, height = height_in, units = "in", res = dpi)
print(p_pcoa)
dev.off()

# --- 4.2. Non-metric Multidimensional Scaling (NMDS) and ANOSIM
nmds_res <- metaMDS(t(data_abund), distance = "bray", k = 2, trymax = 50)

nmds_df <- data.frame(
  Sample = colnames(data_abund),
  Condition = cond,
  NMDS1 = nmds_res$points[,1],
  NMDS2 = nmds_res$points[,2]
)

# ANOSIM for statistical test of differences (associated with NMDS)
anosim_res <- vegan::anosim(dist_bc, grouping = cond, permutations = 999)
write.csv(as.data.frame(anosim_res$class.vec), file.path(out_dir, "ANOSIM_results.csv")) # Save ANOSIM result

# Extract ANOSIM R and p-value for plot annotation
anosim_R <- round(anosim_res$statistic, 3)
anosim_p <- format.pval(anosim_res$signif, digits = 3) # Format p-value for display

# NMDS Plot with ggplot2, including ANOSIM results and connected points
p_nmds <- ggplot(nmds_df, aes(x = NMDS1, y = NMDS2, color = Condition, group = Condition)) + # Added group = Condition
  geom_point(size = 4, aes(shape = Condition)) + # Added shape aesthetic
  stat_ellipse(level = 0.95) + # Add 95% confidence ellipses
  geom_line(aes(group = Condition), alpha = 0.3) + # Conectar puntos dentro de cada grupo
  theme_minimal() +
  labs(title = paste0("NMDS of Gene Abundance Profiles (Stress = ", round(nmds_res$stress, 3), ")"),
       x = "NMDS1",
       y = "NMDS2",
       subtitle = paste0("ANOSIM: R = ", anosim_R, ", p = ", anosim_p)) +
  theme(plot.subtitle = element_text(hjust = 0.5, size = 12)) # Center subtitle

# Save NMDS plot in PNG
png(filename = file.path(out_dir, "NMDS_plot.png"), width = width_in, height = height_in, units = "in", res = dpi)
print(p_nmds)
dev.off()

# Save NMDS plot in TIFF
tiff(filename = file.path(out_dir, "NMDS_plot.tiff"), width = width_in, height = height_in, units = "in", res = dpi)
print(p_nmds)
dev.off()

# ======================
# 5. Variance Partitioning Pie Chart (from PERMANOVA)
# ======================
# Extract R2 and residual
r2_val <- adonis_res$R2[1] # R2 for 'cond'
residual_val <- adonis_res$R2[2] # Residuals

# Create dataframe for pie chart
variance_data <- data.frame(
  Category = c("Explained by Condition", "Residual"),
  Value = c(r2_val, residual_val)
)
variance_data$Percentage <- round(variance_data$Value / sum(variance_data$Value) * 100, 1)

# Pie chart
p_variance_pie <- ggplot(variance_data, aes(x = "", y = Percentage, fill = Category)) +
  geom_bar(width = 1, stat = "identity") +
  coord_polar("y", start = 0) +
  geom_text(aes(label = paste0(Percentage, "%")), position = position_stack(vjust = 0.5), size = 5) +
  labs(title = "Variance Partitioning by Condition (PERMANOVA)",
       fill = "Variance Source") +
  theme_void() + # Minimal theme for pie chart
  theme(plot.title = element_text(hjust = 0.5, size = 16),
        legend.title = element_text(size = 12),
        legend.text = element_text(size = 10))

# Save Pie Chart in PNG
png(filename = file.path(out_dir, "PERMANOVA_variance_pie_chart.png"), width = width_in, height = height_in, units = "in", res = dpi)
print(p_variance_pie)
dev.off()

# Save Pie Chart in TIFF
tiff(filename = file.path(out_dir, "PERMANOVA_variance_pie_chart.tiff"), width = width_in, height = height_in, units = "in", res = dpi)
print(p_variance_pie)
dev.off()


# ======================
# 6. Heatmap of Abundances (ANOVA significant genes)
# ======================

# Max number of genes per condition (can adjust)
N <- 25

# Create a long dataframe for ANOVA directly from data_abund
data_long <- data.frame(
  gene = rep(rownames(data_abund), times = ncol(data_abund)),
  sample = rep(colnames(data_abund), each = nrow(data_abund)),
  abundance = as.vector(as.matrix(data_abund)),
  condition = rep(cond, each = nrow(data_abund))
)

# ANOVA per gene
pvals <- numeric(nrow(data_abund))
names(pvals) <- rownames(data_abund)

for (i in seq_len(nrow(data_abund))) {
  gene_id <- rownames(data_abund)[i]
  sub_data <- subset(data_long, gene == gene_id)
  if (nrow(sub_data) >= 3 && length(unique(sub_data$condition)) > 1) {
    model <- aov(abundance ~ condition, data = sub_data)
    pvals[i] <- summary(model)[[1]][["Pr(>F)"]][1]
  } else {
    pvals[i] <- NA
  }
}

# FDR adjustment
padj <- p.adjust(pvals, method = "fdr")

# Filter significant genes
sig_genes <- names(padj)[which(padj < 0.05 & !is.na(padj))]

# Subset of the original matrix
data_sig <- data_abund[sig_genes, , drop = FALSE]

# If there are significant genes, select the most abundant
if (nrow(data_sig) == 0) {
  warning("No genes with adjusted p-value < 0.05. Heatmap not generated.")
} else {
  # Robust selection: top significant genes with highest mean abundance (global)
  top_genes <- names(sort(rowMeans(data_sig), decreasing = TRUE))[1:min(N, nrow(data_sig))]
  heatmap_data <- data_sig[top_genes, , drop = FALSE]
  
  # Z-score scaling by row
  heatmap_scaled <- t(scale(t(heatmap_data)))
  
  # Annotation by condition
  annotation_col <- data.frame(Condition = cond)
  rownames(annotation_col) <- colnames(data_abund)
  
  # Draw heatmap and save in PNG
  png(filename = file.path(out_dir, "heatmap_ANOVA_top_genes.png"), width = width_in * dpi, height = height_in * dpi, res = dpi)
  pheatmap::pheatmap(
    heatmap_scaled,
    cluster_rows = TRUE,
    cluster_cols = TRUE,
    annotation_col = annotation_col,
    main = "Heatmap of ANOVA Significant Genes" # Added a title
  )
  dev.off()
  
  # Draw heatmap and save in TIFF
  tiff(filename = file.path(out_dir, "heatmap_ANOVA_top_genes.tiff"), width = width_in, height = height_in, units = "in", res = dpi)
  pheatmap::pheatmap(
    heatmap_scaled,
    cluster_rows = TRUE,
    cluster_cols = TRUE,
    annotation_col = annotation_col,
    main = "Heatmap of ANOVA Significant Genes" # Added a title
  )
  dev.off()
}

# ======================
# 7. Venn Diagram of Genes Present in Conditions
# ======================
# Define presence with threshold > 0 (can adjust threshold if needed, e.g. 0.1 CPM)
presence <- data_abund > 0

genes_cond <- list(
  C = rownames(data_abund)[rowSums(presence[, grep("^C_", colnames(data_abund))]) > 0],
  LD = rownames(data_abund)[rowSums(presence[, grep("^LD_", colnames(data_abund))]) > 0],
  PET = rownames(data_abund)[rowSums(presence[, grep("^PET_", colnames(data_abund))]) > 0]
)

# Save Venn diagram in PNG
venn.plot.png <- VennDiagram::venn.diagram(
  x = genes_cond,
  filename = file.path(out_dir, "venn_diagram_genes.png"),
  main = "Genes Present by Condition",
  fill = c("red", "green", "blue"),
  resolution = dpi # Set resolution
)

# Save Venn diagram in TIFF
venn.plot.tiff <- VennDiagram::venn.diagram(
  x = genes_cond,
  filename = file.path(out_dir, "venn_diagram_genes.tiff"),
  main = "Genes Present by Condition",
  fill = c("red", "green", "blue"),
  resolution = dpi # Set resolution
)

# ======================
# 8. Hierarchical Clustering and Dendrogram by Experimental Replicates
# ======================

# Calculate distance between columns (experimental replicates)
dist_replicas <- dist(t(data_abund))

# Hierarchical clustering on replicates
hc_replicas <- hclust(dist_replicas)

# Save dendrogram to PNG file
png(filename = file.path(out_dir, "dendrogram_replicas.png"), width = width_in, height = height_in, units = "in", res = dpi)
plot(hc_replicas, main = "Hierarchical Clustering of Experimental Replicates", xlab = "", sub = "")
dev.off()

# Save dendrogram to TIFF file
tiff(filename = file.path(out_dir, "dendrogram_replicas.tiff"), width = width_in, height = height_in, units = "in", res = dpi)
plot(hc_replicas, main = "Hierarchical Clustering of Experimental Replicates", xlab = "", sub = "")
dev.off()


# ======================
# 9. Gene Co-occurrence Network Analysis
# ======================

# Calculate Spearman correlation between genes (rows)
cor_res <- rcorr(as.matrix(t(data_abund)), type = "spearman") # t() for genes in rows

cor_mat <- cor_res$r         # Matrix of correlation coefficients
p_mat <- cor_res$P           # Matrix of p-values

# Adjust p-values by Benjamini-Hochberg method
p_adj_mat <- matrix(p.adjust(p_mat, method = "BH"), nrow = nrow(p_mat), ncol = ncol(p_mat))
rownames(p_adj_mat) <- rownames(p_mat)
colnames(p_adj_mat) <- colnames(p_mat)

# Filter associations that meet criteria:
# |cor| >= 0.7 and p.adj <= 0.01, excluding self-connections (diagonal)
edges <- which(abs(cor_mat) >= 0.7 & p_adj_mat <= 0.01 & row(cor_mat) != col(cor_mat), arr.ind = TRUE)

# Create edge data frame with source, target, weight, and sign
edge_list <- data.frame(
  from = rownames(cor_mat)[edges[,1]],
  to = colnames(cor_mat)[edges[,2]],
  weight = cor_mat[edges],
  sign = ifelse(cor_mat[edges] > 0, "positive", "negative"),
  stringsAsFactors = FALSE
)

# Create graph without negative weights in layout
g <- graph_from_data_frame(edge_list, directed = FALSE)

# Add absolute weight as attribute for layout (positive)
E(g)$weight_abs <- abs(E(g)$weight)

# Calculate degree (number of connections)
V(g)$degree <- degree(g)

# Define edge colors by sign
E(g)$color <- ifelse(E(g)$sign == "positive", "red", "blue")

# Define node size by degree (scale for better visualization)
V(g)$size <- 5 + 2 * V(g)$degree

# Visualize with Fruchterman-Reingold layout using absolute weights

# Save in PNG
png(filename = file.path(out_dir, "gene_cooccurrence_network.png"), width = width_in, height = height_in, units = "in", res = dpi)
set.seed(123) # For reproducibility of layout
plot(g, layout = layout_with_fr(g, weights = E(g)$weight_abs),
     vertex.label = NA,
     vertex.color = "lightblue",
     main = "Gene Co-occurrence Network (Spearman)")
legend("topright", legend = c("Positive correlation", "Negative correlation"), col = c("red", "blue"), lty = 1, bty = "n")
dev.off()

# Save in TIFF
tiff(filename = file.path(out_dir, "gene_cooccurrence_network.tiff"), width = width_in, height = height_in, units = "in", res = dpi)
set.seed(123) # For reproducibility of layout
plot(g, layout = layout_with_fr(g, weights = E(g)$weight_abs),
     vertex.label = NA,
     vertex.color = "lightblue",
     main = "Gene Co-occurrence Network (Spearman)")
legend("topright", legend = c("Positive correlation", "Negative correlation"), col = c("red", "blue"), lty = 1, bty = "n")
dev.off()

