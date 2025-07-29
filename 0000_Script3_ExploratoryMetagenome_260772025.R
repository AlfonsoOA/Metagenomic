# ====================
# 0. Libraries and Setup
# ====================
packages <- c("tidyverse", "vegan", "ggpubr", "reshape2", "pheatmap",
              "VennDiagram", "ape", "igraph", "psych", "cluster", "dendextend", "Hmisc")

# Install and load packages if not already installed
invisible(lapply(packages, function(x) {
  if (!require(x, character.only = TRUE)) {
    install.packages(x, dependencies = TRUE)
  }
  library(x, character.only = TRUE)
}))

# ===========================
# 1. Load and Prepare Data
# ===========================
base_path <- "C:/Users/Alfonso/Desktop/AlfonsoOA_MSI/1Investigacion/5BiodegradacionPlasticos/11MetagenomaTenebrio/0_Metagenome_Figure2"
file <- file.path(base_path, "1_EC_relative_abundance_table_tenebrio_molitor.csv")
out_dir <- file.path(base_path, "Metagenome_Results")

# Ensure output directory exists
if (!dir.exists(out_dir)) {
  dir.create(out_dir, recursive = TRUE)
}

# Load data without factors to avoid issues
data_raw <- read.csv(file, sep = "\t", stringsAsFactors = FALSE)

# Set first column as row names (assuming 'gene' column is the first)
rownames(data_raw) <- data_raw$gene

# Select abundance columns (columns 5 to 16 = 3 conditions * 4 replicates)
abund_cols <- colnames(data_raw)[5:16]

# Extract abundance matrix: rows = genes, columns = samples
data_abund <- data_raw[, abund_cols]

# Convert to numeric matrix for analysis
data_abund <- as.matrix(data_abund)
mode(data_abund) <- "numeric"

# Define GLOBAL condition factor based on data_abund sample names
# This 'cond_global' will be used for sections 2, 4, 5, 6, 7, 8
cond_global <- factor(gsub("_\\dK", "", abund_cols), levels = c("C", "LD", "PET"))

# Configuration for saving high-quality images
dpi <- 300 # Kept at 300 DPI
width_in <- 7
height_in <- 5

# ==================
# 2. Alpha Diversity
# ==================
alpha_div <- data.frame(
  sample = colnames(data_abund),
  condition = cond_global, # Use global condition
  Shannon = vegan::diversity(t(data_abund), index = "shannon"),
  Simpson = vegan::diversity(t(data_abund), index = "simpson")
)

write.csv(alpha_div, file.path(out_dir, "alpha_diversity.csv"), row.names = FALSE)

# Function to create and save alpha diversity plot
plot_alpha_diversity <- function(data, index_name, out_dir, width, height, res) {
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
  tiff(tiff_filename, width = width, height = height, units = "in", res = res)
  print(p)
  dev.off()
  
  # Save in PNG
  png_filename <- file.path(out_dir, paste0("alpha_diversity_", index_name, ".png"))
  png(png_filename, width = width, height = height, units = "in", res = res)
  print(p)
  dev.off()
}

# Call function with global dpi
plot_alpha_diversity(alpha_div, "Shannon", out_dir, width_in, height_in, dpi)
plot_alpha_diversity(alpha_div, "Simpson", out_dir, width_in, height_in, dpi)

# =====================
# 3. Sample Rarefaction (from Salmon counts)
# =====================

# Directory where your .sf files are located
salmon_dir <- file.path(base_path, "Tenebrio_molitor_Salmon") 

# Get list of _quant.sf files
quant_files <- list.files(salmon_dir, pattern = "_quant.sf$", full.names = TRUE)

# Filter out the "_sub_quant.sf" files
quant_files <- quant_files[!grepl("_sub_quant.sf$", quant_files)]

# Filter out empty or unreadable files if any
valid_quant_files <- c()
for (f in quant_files) {
  if (file.exists(f) && file.size(f) > 0) {
    valid_quant_files <- c(valid_quant_files, f)
  } else {
    warning(paste("Skipping empty or non-existent Salmon file:", f))
  }
}
quant_files <- valid_quant_files

if (length(quant_files) == 0) {
  stop("No valid Salmon quant.sf files found (excluding _sub_quant.sf). Cannot perform rarefaction.")
}

# Extract sample names from file names (e.g., C_1K_MaxBin2Refined-0.002_quant.sf)
sample_names_salmon <- basename(quant_files)

# Initialize a list to store counts for each sample
all_counts_list <- list()

# Loop to read each .sf file
for (i in seq_along(quant_files)) {
  file_path <- quant_files[i]
  sample_name <- sample_names_salmon[i]
  
  # Read only the Name and NumReads columns
  df_salmon <- read.table(file_path, header = TRUE, sep = "\t", quote = "", stringsAsFactors = FALSE)
  
  # Ensure 'NumReads' are integers (they are expected counts from Salmon)
  # round() is good practice for rarefaction functions.
  all_counts_list[[sample_name]] <- setNames(round(df_salmon$NumReads), df_salmon$Name)
}

# Consolidate all counts into a single matrix
# First, find all unique gene/transcript names
all_genes <- unique(unlist(lapply(all_counts_list, names)))

# Create the final count matrix, filling with 0 if a gene is not present in a sample
raw_counts_matrix <- matrix(0, nrow = length(all_genes), ncol = length(sample_names_salmon),
                            dimnames = list(all_genes, sample_names_salmon))

for (sample_name in sample_names_salmon) {
  current_counts <- all_counts_list[[sample_name]]
  # Ensure genes from current_counts are in all_genes before assigning
  valid_genes_in_current <- names(current_counts)[names(current_counts) %in% all_genes]
  if(length(valid_genes_in_current) > 0) {
    raw_counts_matrix[valid_genes_in_current, sample_name] <- current_counts[valid_genes_in_current]
  }
}

# Transpose the matrix so rows are samples and columns are genes
# (format expected by vegan's rarecurve function)
data_for_rare <- t(raw_counts_matrix)

# Remove samples with zero total reads from data_for_rare, as they cause issues with rarefaction
initial_samples_count <- nrow(data_for_rare)
data_for_rare <- data_for_rare[rowSums(data_for_rare) > 0, , drop = FALSE]
if (nrow(data_for_rare) < initial_samples_count) {
  warning(paste(initial_samples_count - nrow(data_for_rare), "samples removed from rarefaction data because they had zero total reads."))
}

# Ensure all counts are non-negative (should be, after round and matrix initialization)
stopifnot(all(data_for_rare >= 0))

# Define condition factor for the Salmon data used in rarefaction.
# Use a more robust pattern for sample names like "C_1K_MaxBin2Refined-0.002_quant.sf"
# This pattern extracts 'C', 'LD', or 'PET' from the beginning of the string.
# It should now correctly handle the *filtered* list of names.
cond_salmon <- factor(gsub("^(C|LD|PET)_.*", "\\1", rownames(data_for_rare)), levels = c("C", "LD", "PET"))

# Check for NAs in cond_salmon after its creation. If NAs, it means the gsub pattern failed for some samples.
if (any(is.na(cond_salmon))) {
  warning("NA values found in 'cond_salmon'. This indicates that some sample names did not match the expected pattern 'C_.*', 'LD_.*', 'PET_.*'. Please check your Salmon sample names. Filtering samples with NA condition to avoid errors.")
  # Filter out samples with NA condition to avoid errors in plotting.
  data_for_rare <- data_for_rare[!is.na(cond_salmon), , drop = FALSE]
  cond_salmon <- cond_salmon[!is.na(cond_salmon)]
}

# Colors by condition
cols_cond <- c(C = "blue", LD = "red", PET = "green")

# Minimum total reads per sample for rarefaction saturation point
raremax <- if(nrow(data_for_rare) > 0) min(rowSums(data_for_rare)) else 0

# Debugging check: if raremax is 0, rarecurve will likely be empty.
if (raremax == 0 || nrow(data_for_rare) == 0) {
  warning("Minimum number of reads for rarefaction (raremax) is 0 or no samples remain after filtering. Rarefaction curves were not generated.")
} else {
  # Estimate max number of genes that could be detected for ylim
  max_genes_detected_estimate <- max(apply(data_for_rare, 1, function(x) sum(x > 0))) * 1.1 
  if (is.infinite(max_genes_detected_estimate) || max_genes_detected_estimate == 0) max_genes_detected_estimate <- 100 # Fallback if data is too sparse
  
  # Save rarefaction curve in TIFF
  tryCatch({
    tiff(filename = file.path(out_dir, "Rarefaction_curve.tiff"), width = width_in, height = height_in, units = "in", res = dpi)
    rarecurve(data_for_rare, step = 100, sample = raremax, label = FALSE,
              col = cols_cond[cond_salmon], lwd = 2,
              xlab = "Reads per sample", ylab = "Genes detected", main = "Rarefaction Curves",
              ylim = c(0, max_genes_detected_estimate))
    legend("bottomright", legend = levels(cond_salmon), col = cols_cond, lwd = 2, bty = "n")
    dev.off()
  }, error = function(e) {
    message("Error generating Rarefaction_curve.tiff: ", e$message)
    if (!is.null(dev.list())) dev.off()
  })
  
  # Save rarefaction curve in PNG
  tryCatch({
    png(filename = file.path(out_dir, "Rarefaction_curve.png"), width = width_in, height = height_in, units = "in", res = dpi)
    rarecurve(data_for_rare, step = 100, sample = raremax, label = FALSE,
              col = cols_cond[cond_salmon], lwd = 2,
              xlab = "Reads per sample", ylab = "Genes detected", main = "Rarefaction Curves",
              ylim = c(0, max_genes_detected_estimate))
    legend("bottomright", legend = levels(cond_salmon), col = cols_cond, lwd = 2, bty = "n")
    dev.off()
  }, error = function(e) {
    message("Error generating Rarefaction_curve.png: ", e$message)
    if (!is.null(dev.list())) dev.off()
  })
}

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
  Condition = cond_global, # Use global condition
  Axis1 = pcoa_res$vectors[,1],
  Axis2 = pcoa_res$vectors[,2]
)

# Save PCoA table
write.csv(pcoa_df, file.path(out_dir, "PCoA_coordinates.csv"), row.names = FALSE)

# PERMANOVA for statistical test of differences (associated with PCoA)
# This uses cond_global directly from Section 1, which has been checked for NAs.
adonis_res <- vegan::adonis2(dist_bc ~ cond_global, permutations = 999) 
write.csv(as.data.frame(adonis_res), file.path(out_dir, "PERMANOVA_results.csv"))

# Extract PERMANOVA R2 and p-value for plot annotation
permanova_r2 <- round(adonis_res$R2[1], 3)
permanova_p <- format.pval(adonis_res$`Pr(>F)`[1], digits = 3) # Format p-value for display

# PCoA Plot with ggplot2, including PERMANOVA results and connected points
p_pcoa <- ggplot(pcoa_df, aes(x = Axis1, y = Axis2, color = Condition, group = Condition)) + 
  geom_point(size = 4, aes(shape = Condition)) + 
  stat_ellipse(level = 0.95) + 
  geom_line(aes(group = Condition), alpha = 0.3) + 
  theme_minimal() +
  labs(title = "PCoA of Gene Abundance Profiles",
       x = paste0("PCoA1 (", round(pcoa_res$values$Relative_eig[1]*100, 2), "%)"),
       y = paste0("PCoA2 (", round(pcoa_res$values$Relative_eig[2]*100, 2), "%)"),
       subtitle = paste0("PERMANOVA: R² = ", permanova_r2, ", p = ", permanova_p)) +
  theme(plot.subtitle = element_text(hjust = 0.5, size = 12)) 

# Save PCoA plot in PNG
tryCatch({
  png(filename = file.path(out_dir, "PCoA_plot.png"), width = width_in, height = height_in, units = "in", res = dpi)
  print(p_pcoa)
  dev.off()
}, error = function(e) {
  message("Error generating PCoA_plot.png: ", e$message)
  if (!is.null(dev.list())) dev.off()
})

# Save PCoA plot in TIFF
tryCatch({
  tiff(filename = file.path(out_dir, "PCoA_plot.tiff"), width = width_in, height = height_in, units = "in", res = dpi)
  print(p_pcoa)
  dev.off()
}, error = function(e) {
  message("Error generating PCoA_plot.tiff: ", e$message)
  if (!is.null(dev.list())) dev.off()
})

# --- 4.2. Non-metric Multidimensional Scaling (NMDS) and ANOSIM
nmds_res <- metaMDS(t(data_abund), distance = "bray", k = 2, trymax = 50)

nmds_df <- data.frame(
  Sample = colnames(data_abund),
  Condition = cond_global, # Use global condition
  NMDS1 = nmds_res$points[,1],
  NMDS2 = nmds_res$points[,2]
)

# ANOSIM for statistical test of differences (associated with NMDS)
# This uses cond_global directly from Section 1, which has been checked for NAs.
anosim_res <- vegan::anosim(dist_bc, grouping = cond_global, permutations = 999) 
write.csv(as.data.frame(anosim_res$class.vec), file.path(out_dir, "ANOSIM_results.csv")) 

# Extract ANOSIM R and p-value for plot annotation
anosim_R <- round(anosim_res$statistic, 3)
anosim_p <- format.pval(anosim_res$signif, digits = 3) 

# NMDS Plot with ggplot2, including ANOSIM results and connected points
p_nmds <- ggplot(nmds_df, aes(x = NMDS1, y = NMDS2, color = Condition, group = Condition)) + 
  geom_point(size = 4, aes(shape = Condition)) + 
  stat_ellipse(level = 0.95) + 
  geom_line(aes(group = Condition), alpha = 0.3) + 
  theme_minimal() +
  labs(title = paste0("NMDS of Gene Abundance Profiles (Stress = ", round(nmds_res$stress, 3), ")"),
       x = "NMDS1",
       y = "NMDS2",
       subtitle = paste0("ANOSIM: R = ", anosim_R, ", p = ", anosim_p)) +
  theme(plot.subtitle = element_text(hjust = 0.5, size = 12)) 

# Save NMDS plot in PNG
tryCatch({
  png(filename = file.path(out_dir, "NMDS_plot.png"), width = width_in, height = height_in, units = "in", res = dpi)
  print(p_nmds)
  dev.off()
}, error = function(e) {
  message("Error generating NMDS_plot.png: ", e$message)
  if (!is.null(dev.list())) dev.off()
})

# Save NMDS plot in TIFF
tryCatch({
  tiff(filename = file.path(out_dir, "NMDS_plot.tiff"), width = width_in, height = height_in, units = "in", res = dpi)
  print(p_nmds)
  dev.off()
}, error = function(e) {
  message("Error generating NMDS_plot.tiff: ", e$message)
  if (!is.null(dev.list())) dev.off()
})

# ======================
# 5. Variance Partitioning Pie Chart (from PERMANOVA)
# ======================
# Extract R2 and residual
r2_val <- adonis_res$R2[1] # R2 for 'cond_global'
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
  theme_void() + 
  theme(plot.title = element_text(hjust = 0.5, size = 16),
        legend.title = element_text(size = 12),
        legend.text = element_text(size = 10))

# Save Pie Chart in PNG
tryCatch({
  png(filename = file.path(out_dir, "PERMANOVA_variance_pie_chart.png"), width = width_in, height = height_in, units = "in", res = dpi)
  print(p_variance_pie)
  dev.off()
}, error = function(e) {
  message("Error generating PERMANOVA_variance_pie_chart.png: ", e$message)
  if (!is.null(dev.list())) dev.off()
})

# Save Pie Chart in TIFF
tryCatch({
  tiff(filename = file.path(out_dir, "PERMANOVA_variance_pie_chart.tiff"), width = width_in, height = height_in, units = "in", res = dpi)
  print(p_variance_pie)
  dev.off()
}, error = function(e) {
  message("Error generating PERMANOVA_variance_pie_chart.tiff: ", e$message)
  if (!is.null(dev.list())) dev.off()
})

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
  condition = rep(cond_global, each = nrow(data_abund)) # Use global condition
)

# ANOVA per gene
pvals <- numeric(nrow(data_abund))
names(pvals) <- rownames(data_abund)

for (i in seq_len(nrow(data_abund))) {
  gene_id <- rownames(data_abund)[i]
  sub_data <- subset(data_long, gene == gene_id)
  
  # Ensure sufficient data for ANOVA (at least 2 observations per group for 2 groups, 3 for 3 groups)
  # And at least one unique condition.
  if (nrow(sub_data) >= length(unique(sub_data$condition)) && length(unique(sub_data$condition)) > 1) {
    group_variances <- tapply(sub_data$abundance, sub_data$condition, var)
    
    if (any(is.na(group_variances))) {
      pvals[i] <- NA
    } else if (all(group_variances > 0)) {
      model <- aov(abundance ~ condition, data = sub_data)
      pvals[i] <- summary(model)[[1]][["Pr(>F)"]][1]
    } else {
      pvals[i] <- NA
    }
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
  mean_abundances_sig <- rowMeans(data_sig)
  
  # Robust selection: top significant genes with highest mean abundance (global)
  top_genes <- names(sort(mean_abundances_sig, decreasing = TRUE))[1:min(N, nrow(data_sig))]
  heatmap_data <- data_sig[top_genes, , drop = FALSE]
  
  # Check if there are enough genes and samples for clustering for pheatmap
  if (nrow(heatmap_data) < 2 || ncol(heatmap_data) < 2) {
    warning("Not enough genes or samples (less than 2) selected for heatmap clustering. Heatmap not generated.")
  } else {
    # Z-score scaling by row
    row_variances <- apply(heatmap_data, 1, var)
    if (any(row_variances == 0)) {
      warning("Some selected genes have zero variance across samples. Setting their scaled values to 0.")
      heatmap_scaled <- t(apply(heatmap_data, 1, function(x) {
        if(var(x) == 0) rep(0, length(x)) else scale(x)
      }))
    } else {
      heatmap_scaled <- t(scale(t(heatmap_data)))
    }
    
    # Force garbage collection before plotting
    gc() 
    
    # Annotation by condition
    annotation_col <- data.frame(Condition = cond_global) 
    rownames(annotation_col) <- colnames(data_abund)
    
    # Try/catch block to ensure dev.off() is called even if an error occurs
    tryCatch({
      # Save in PNG
      png(filename = file.path(out_dir, "heatmap_ANOVA_top_genes.png"), 
          width = width_in, height = height_in, units = "in", res = dpi)
      # Ensure all NA/Inf/NaN are handled for pheatmap
      heatmap_scaled[is.na(heatmap_scaled)] <- 0
      heatmap_scaled[is.infinite(heatmap_scaled)] <- 0
      pheatmap::pheatmap(
        heatmap_scaled,
        cluster_rows = TRUE,
        cluster_cols = TRUE,
        annotation_col = annotation_col,
        main = "Heatmap of ANOVA Significant Genes"
      )
      dev.off()
      
      # Save in TIFF
      tiff(filename = file.path(out_dir, "heatmap_ANOVA_top_genes.tiff"), 
           width = width_in, height = height_in, units = "in", res = dpi)
      # Ensure all NA/Inf/NaN are handled for pheatmap
      heatmap_scaled[is.na(heatmap_scaled)] <- 0
      heatmap_scaled[is.infinite(heatmap_scaled)] <- 0
      pheatmap::pheatmap(
        heatmap_scaled,
        cluster_rows = TRUE,
        cluster_cols = TRUE,
        annotation_col = annotation_col,
        main = "Heatmap of ANOVA Significant Genes"
      )
      dev.off()
      
    }, error = function(e) {
      message("Error generating heatmap: ", e$message)
      if (!is.null(dev.list())) {
        dev.off()
      }
    })
  } 
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
tryCatch({
  venn.plot.png <- VennDiagram::venn.diagram(
    x = genes_cond,
    filename = file.path(out_dir, "venn_diagram_genes.png"),
    main = "Genes Present by Condition",
    fill = c("red", "green", "blue"),
    resolution = dpi 
  )
}, error = function(e) {
  message("Error generating venn_diagram_genes.png: ", e$message)
})


# Save Venn diagram in TIFF
tryCatch({
  venn.plot.tiff <- VennDiagram::venn.diagram(
    x = genes_cond,
    filename = file.path(out_dir, "venn_diagram_genes.tiff"),
    main = "Genes Present by Condition",
    fill = c("red", "green", "blue"),
    resolution = dpi 
  )
}, error = function(e) {
  message("Error generating venn_diagram_genes.tiff: ", e$message)
})

# ======================
# 8. Hierarchical Clustering and Dendrogram by Experimental Replicates
# ======================

# Calculate distance between columns (experimental replicates)
dist_replicas <- dist(t(data_abund))

# Hierarchical clustering on replicates
hc_replicas <- hclust(dist_replicas)

# Save dendrogram to PNG file
tryCatch({
  png(filename = file.path(out_dir, "dendrogram_replicas.png"), width = width_in, height = height_in, units = "in", res = dpi)
  plot(hc_replicas, main = "Hierarchical Clustering of Experimental Replicates", xlab = "", sub = "")
  dev.off()
}, error = function(e) {
  message("Error generating dendrogram_replicas.png: ", e$message)
  if (!is.null(dev.list())) dev.off()
})

# Save dendrogram to TIFF file
tryCatch({
  tiff(filename = file.path(out_dir, "dendrogram_replicas.tiff"), width = width_in, height = height_in, units = "in", res = dpi)
  plot(hc_replicas, main = "Hierarchical Clustering of Experimental Replicates", xlab = "", sub = "")
  dev.off()
}, error = function(e) {
  message("Error generating dendrogram_replicas.tiff: ", e$message)
  if (!is.null(dev.list())) dev.off()
})

# ======================
# 9. Gene Co-occurrence Network Analysis
# ======================

# Calculate Spearman correlation between genes (rows)
cor_res <- rcorr(as.matrix(t(data_abund)), type = "spearman") 

cor_mat <- cor_res$r         
p_mat <- cor_res$P           

# Adjust p-values by Benjamini-Hochberg method
p_adj_mat <- matrix(p.adjust(p_mat, method = "BH"), nrow = nrow(p_mat), ncol = ncol(p_mat))
rownames(p_adj_mat) <- rownames(p_mat)
colnames(p_adj_mat) <- colnames(p_mat)

# Filter associations that meet criteria:
# Adjust these based on biological significance and network density
threshold_corr <- 0.6 
threshold_padj <- 0.05 

edges <- which(abs(cor_mat) >= threshold_corr & p_adj_mat <= threshold_padj & row(cor_mat) != col(cor_mat), arr.ind = TRUE)

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

# Only attempt to plot if graph has vertices (nodes)
if (vcount(g) > 0) {
  # Save in PNG
  tryCatch({
    png(filename = file.path(out_dir, "gene_cooccurrence_network.png"), width = width_in, height = height_in, units = "in", res = dpi)
    set.seed(123) # For reproducibility of layout
    plot(g, layout = layout_with_fr(g, weights = E(g)$weight_abs),
         vertex.label = NA,
         vertex.color = "lightblue",
         main = "Gene Co-occurrence Network (Spearman)")
    legend("topright", legend = c("Positive correlation", "Negative correlation"), col = c("red", "blue"), lty = 1, bty = "n")
    dev.off()
  }, error = function(e) {
    message("Error generating gene_cooccurrence_network.png: ", e$message)
    if (!is.null(dev.list())) dev.off()
  })
  
  # Save in TIFF
  tryCatch({
    tiff(filename = file.path(out_dir, "gene_cooccurrence_network.tiff"), width = width_in, height = height_in, units = "in", res = dpi)
    set.seed(123) # For reproducibility of layout
    plot(g, layout = layout_with_fr(g, weights = E(g)$weight_abs),
         vertex.label = NA,
         vertex.color = "lightblue",
         main = "Gene Co-occurrence Network (Spearman)")
    legend("topright", legend = c("Positive correlation", "Negative correlation"), col = c("red", "blue"), lty = 1, bty = "n")
    dev.off()
  }, error = function(e) {
    message("Error generating gene_cooccurrence_network.tiff: ", e$message)
    if (!is.null(dev.list())) dev.off()
  })
} else {
  warning("No nodes found in gene co-occurrence network (no significant correlations). Network plot not generated.")
}

