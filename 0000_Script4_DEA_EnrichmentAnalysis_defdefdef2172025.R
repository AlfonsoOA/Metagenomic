# Load necessary libraries
# install.packages(c("tidyverse", "readr", "dplyr", "ggplot2", "tidyr", "scales", "forcats", "limma", "clusterProfiler", "DOSE", "enrichplot", "ggrepel", "stringr", "cowplot"))

library(tidyverse)
library(readr)
library(dplyr)
library(ggplot2)
library(tidyr)
library(scales) # For scales::label_number
library(forcats) # For fct_reorder
library(limma)
library(clusterProfiler)
library(DOSE) # Often loaded with clusterProfiler
library(enrichplot) # For dotplot
library(ggrepel) # For better label placement
library(stringr) # For string manipulation
library(cowplot) # For plot_grid

# ==============================================================================
# 1. Configuration & Global Parameters ⚙️
#    (Adjust these parameters as needed)
# ==============================================================================

# Base path for input and output files
BASE_PATH <- "C:/Users/Alfonso/Desktop/AlfonsoOA_MSI/1Investigacion/5BiodegradacionPlasticos/11MetagenomaTenebrio/00_AnotacionInterProScan+GhostKOALA"

# Input file paths
FILE_PATHS <- list(
  ec_abundance = file.path(BASE_PATH, "1_EC_relative_abundance_table_tenebrio_molitor.csv"),
  omicbox_annotation = file.path(BASE_PATH, "OmicBox_annotation.csv"),
  interproscan_full = file.path(BASE_PATH, "00_MetaBAT2Refined-0.2_interproscan.tsv"),
  ko_annotations = file.path(BASE_PATH, "user_ko.txt")
)

# Output directory for all results
OUTPUT_RESULTS_DIR <- file.path(BASE_PATH, "AnnotationPGresults_v3")

# Plot dimensions and resolution
PLOT_WIDTH_CM <- 25
PLOT_HEIGHT_CM <- 20 # Adjusted for better fit
DPI <- 600

# Expected column names and patterns
PROTEIN_ID_COL_NAME_MAIN <- "gene" # Column name for protein accession in main_data_df
CPM_COLS_PATTERN <- "^(C|LD|PET)_\\dK$" # Pattern for CPM columns (e.g., C_1K, LD_2K)

# Delimiters for input files (CRITICAL: Adjust if your files use commas or spaces)
DELIMITERS <- list(
  ec_abundance = ";",
  omicbox_annotation = ";",
  interproscan_full = "\t", # InterProScan is typically tab-separated
  ko_annotations = "\t" # user_ko.txt is typically tab-separated
)

# Threshold for pre-DEA Fold Change categorization (log2FC)
PRE_DEA_LOG2FC_THRESHOLD <- 1 # E.g., 1 means 2-fold change

# Pseudocount for log2 transformation in limma and plot visualization
PSEUDOCOUNT <- 1e-6 # A very small number to avoid log(0) errors

# ==============================================================================
# 2. Helper Functions for Annotation Analysis 🛠️
# ==============================================================================

#' Plots abundance of InterProScan methods or other categories by condition.
#'
#' @param data A dataframe with Method and CPM_Value columns.
#' @param output_dir Directory to save plots.
#' @param file_suffix Suffix for output file names.
#' @param width_cm Plot width in cm.
#' @param height_cm Plot height in cm.
#' @param dpi Plot resolution.
plot_abundance_by_method <- function(data, output_dir, file_suffix, width_cm, height_cm, dpi) {
  df_long <- data %>%
    pivot_longer(
      cols = starts_with(c("C_", "LD_", "PET_")),
      names_to = "Sample",
      values_to = "CPM_Value"
    ) %>%
    mutate(Condition = sub("_[0-9]+K", "", Sample))
  
  summary_df <- df_long %>%
    group_by(Method, Condition) %>%
    summarise(Total_CPM = sum(CPM_Value, na.rm = TRUE), .groups = 'drop') %>%
    filter(Total_CPM > 0)
  
  if (nrow(summary_df) == 0) {
    cat(paste0("No data to plot for abundance by method (", file_suffix, "). Skipping plot.\n"))
    return(NULL)
  }
  
  summary_df$Method <- fct_reorder(summary_df$Method, summary_df$Total_CPM, .fun = sum)
  
  p <- ggplot(summary_df, aes(x = Method, y = Total_CPM + PSEUDOCOUNT, fill = Condition)) + # Added pseudocount for plotting
    geom_bar(stat = "identity", position = position_dodge(width = 0.8), width = 0.7) +
    labs(
      title = paste0("Abundance of InterProScan Methods (", file_suffix, ")"),
      x = "InterProScan Method",
      y = "Total CPM (log10)" # Updated label to reflect log scale
    ) +
    scale_y_log10(labels = scales::label_number(scale_cut = scales::cut_short_scale())) +
    scale_fill_brewer(palette = "Set1") +
    theme_minimal(base_size = 10) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 12),
      axis.title = element_text(size = 9),
      axis.text.x = element_text(angle = 45, hjust = 1, size = 7),
      axis.text.y = element_text(size = 8),
      legend.position = "right",
      legend.title = element_text(size = 9, face = "bold"),
      panel.grid.major.x = element_blank(),
      panel.grid.minor.x = element_blank()
    )
  
  output_tiff <- file.path(output_dir, paste0("1_Abundancia_por_metodo_", file_suffix, ".tiff"))
  output_png <- file.path(output_dir, paste0("1_Abundancia_por_metodo_", file_suffix, ".png"))
  ggsave(output_tiff, plot = p, width = width_cm, height = height_cm, units = "cm", dpi = dpi, compression = "lzw")
  ggsave(output_png, plot = p, width = width_cm, height = height_cm, units = "cm", dpi = dpi)
  cat(paste0("Abundance plot saved to: ", output_tiff, " and ", output_png, "\n"))
  return(p)
}

#' Plots cumulative unique annotation counts (rarefaction-like).
#'
#' @param data A dataframe with Protein_Accession and the annotation column.
#' @param annotation_col The name of the column containing annotations (e.g., "Annotation.GO.ID", "Analysis_Method").
#' @param annotation_type_name A descriptive name for the annotation type (e.g., "GO Terms", "InterPro Methods").
#' @param output_dir Directory to save plots.
#' @param file_suffix Suffix for output file names.
#' @param width_cm Plot width in cm.
#' @param height_cm Plot height in cm.
#' @param dpi Plot resolution.
plot_cumulative_annotations <- function(data, annotation_col, annotation_type_name, output_dir, file_suffix, width_cm, height_cm, dpi) {
  
  # Filter out NA/empty annotations and ensure distinct protein-annotation pairs
  protein_annotations <- data %>%
    dplyr::select(Protein_Accession, !!sym(annotation_col)) %>%
    filter(!is.na(!!sym(annotation_col)) & !!sym(annotation_col) != "") %>% # Ensure non-NA and non-empty
    distinct() %>% # Important: Ensures each protein-annotation combination is counted once
    arrange(Protein_Accession) # Order by protein accession for consistent cumulative count
  
  if (nrow(protein_annotations) == 0) {
    cat(paste0("No data to plot for cumulative unique annotations (", annotation_type_name, " - ", file_suffix, "). Skipping plot.\n"))
    return(NULL)
  }
  
  # Calculate cumulative unique annotations
  unique_annotations_cumulative <- tibble(
    Protein_Count = 1:nrow(protein_annotations),
    Unique_Annotations_Seen = sapply(1:nrow(protein_annotations), function(i) {
      length(unique(protein_annotations[[annotation_col]][1:i]))
    })
  )
  
  p <- ggplot(unique_annotations_cumulative, aes(x = Protein_Count, y = Unique_Annotations_Seen)) +
    geom_line(color = "darkblue", linewidth = 0.8) +
    geom_point(color = "darkblue", size = 1.5) +
    labs(
      title = paste0("Cumulative Unique ", annotation_type_name, " (", file_suffix, ")"),
      x = "Number of Proteins Processed (Ordered by Accession)",
      y = paste0("Cumulative Unique ", annotation_type_name)
    ) +
    theme_minimal(base_size = 10) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 12),
      axis.title = element_text(size = 9),
      axis.text = element_text(size = 8)
    )
  
  # Simplified output file names for rarefaction curves
  output_tiff <- file.path(output_dir, paste0("2_Rarefaction_", gsub(" ", "_", file_suffix), ".tiff"))
  output_png <- file.path(output_dir, paste0("2_Rarefaction_", gsub(" ", "_", file_suffix), ".png"))
  ggsave(output_tiff, plot = p, width = width_cm, height = height_cm, units = "cm", dpi = dpi, compression = "lzw")
  ggsave(output_png, plot = p, width = width_cm, height = height_cm, units = "cm", dpi = dpi)
  cat(paste0("Cumulative unique methods plot saved to: ", output_tiff, " and ", output_png, "\n"))
  return(p)
}


#' Plots accumulated CPMs for a given annotation type by condition.
#'
#' @param data_df The dataframe containing CPMs and annotations.
#' @param cpm_cols Names of CPM columns.
#' @param annotation_col Name of the annotation column to plot.
#' @param annotation_type_name Display name for the annotation type.
#' @param output_dir Directory to save plots.
#' @param dpi Plot resolution.
#' @param plot_type "global_top15" for overall top 15, or a specific condition like "C", "LD", "PET" for top 15 of that condition.
#' @param granularity_label An additional label for granular plots (e.g., "Biological Process", "Pfam").
plot_accumulated_cpm <- function(data_df, cpm_cols, annotation_col, annotation_type_name, output_dir, dpi, plot_type = "global_top15", granularity_label = NULL) {
  
  # Ensure annotation_col exists in data_df before proceeding
  if (!(annotation_col %in% colnames(data_df))) {
    cat(paste0("  Error: Column '", annotation_col, "' not found in data for accumulated CPM plot (", 
               ifelse(!is.null(granularity_label), paste0(" (", granularity_label, ")"), ""), 
               "). Skipping plot.\n"))
    return(NULL)
  }
  
  df_filtered <- data_df %>%
    filter(!is.na(!!sym(annotation_col)) & !!sym(annotation_col) != "") # Ensure non-NA and non-empty
  
  if (nrow(df_filtered) == 0) {
    cat(paste0("  No valid annotations for ", annotation_type_name, 
               ifelse(!is.null(granularity_label), paste0(" (", granularity_label, ")"), ""), 
               ". Skipping accumulated CPM plot for ", plot_type, ".\n"))
    return(NULL)
  }
  
  cpm_long <- df_filtered %>%
    dplyr::select(Protein_Accession, !!sym(annotation_col), all_of(cpm_cols)) %>%
    pivot_longer(
      cols = all_of(cpm_cols),
      names_to = "Sample",
      values_to = "CPM_Value"
    ) %>%
    mutate(Condition = sub("_[0-9]+K", "", Sample)) %>%
    group_by(Condition, Sample, Annotation = !!sym(annotation_col)) %>%
    summarise(Accumulated_CPM = sum(CPM_Value, na.rm = TRUE), .groups = 'drop')
  
  # Filter out annotations with 0 total CPMs across all conditions for plotting
  cpm_long <- cpm_long %>%
    group_by(Annotation) %>%
    filter(sum(Accumulated_CPM) > 0) %>%
    ungroup() %>%
    filter(Accumulated_CPM > 0) # Filter individual 0 CPMs for log scale visualization
  
  if (nrow(cpm_long) == 0) {
    cat(paste0("  No non-zero accumulated CPMs for ", annotation_type_name, 
               ifelse(!is.null(granularity_label), paste0(" (", granularity_label, ")"), ""), 
               ". Skipping plot for ", plot_type, ".\n"))
    return(NULL)
  }
  
  plot_title_base <- paste0("Top 15 Accumulated CPMs for ", annotation_type_name)
  output_file_name_suffix <- ""
  selected_annotations <- character(0)
  
  if (plot_type == "global_top15") {
    plot_title <- paste0(plot_title_base, " (Overall)")
    output_file_name_suffix <- "Overall"
    
    selected_annotations <- cpm_long %>%
      group_by(Annotation) %>%
      summarise(Total_Accumulated_CPM = sum(Accumulated_CPM, na.rm = TRUE), .groups = 'drop') %>%
      arrange(desc(Total_Accumulated_CPM)) %>%
      head(15) %>%
      pull(Annotation)
    
  } else if (plot_type %in% c("C", "LD", "PET")) {
    plot_title <- paste0(plot_title_base, " in Condition '", plot_type, "'")
    output_file_name_suffix <- paste0("Top15_in_", plot_type)
    
    selected_annotations <- cpm_long %>%
      filter(Condition == plot_type) %>%
      group_by(Annotation) %>%
      summarise(Total_Accumulated_CPM = sum(Accumulated_CPM, na.rm = TRUE), .groups = 'drop') %>%
      arrange(desc(Total_Accumulated_CPM)) %>%
      head(15) %>%
      pull(Annotation)
    
    if (length(selected_annotations) == 0) {
      cat(paste0("  No top 15 annotations found for '", plot_type, "' condition for ", annotation_type_name, 
                 ifelse(!is.null(granularity_label), paste0(" (", granularity_label, ")"), ""), ". Skipping plot.\n"))
      return(NULL)
    }
    
  } else {
    stop("Invalid plot_type provided. Must be 'global_top15', 'C', 'LD', or 'PET'.")
  }
  
  cpm_long_filtered <- cpm_long %>%
    filter(Annotation %in% selected_annotations)
  
  if (nrow(cpm_long_filtered) == 0) {
    cat(paste0("  No annotations remain after selecting top 15 for ", plot_type, " for ", annotation_type_name, 
               ifelse(!is.null(granularity_label), paste0(" (", granularity_label, ")"), ""), ". Skipping plot.\n"))
    return(NULL)
  }
  
  if (plot_type %in% c("C", "LD", "PET")) {
    order_annotations <- cpm_long_filtered %>%
      filter(Condition == plot_type) %>%
      group_by(Annotation) %>%
      summarise(Mean_CPM_In_Focus_Condition = mean(Accumulated_CPM, na.rm = TRUE), .groups = 'drop') %>%
      arrange(desc(Mean_CPM_In_Focus_Condition)) %>%
      pull(Annotation)
  } else { # For global_top15
    order_annotations <- cpm_long_filtered %>%
      group_by(Annotation) %>%
      summarise(Total_CPM_Overall = sum(Accumulated_CPM, na.rm = TRUE), .groups = 'drop') %>%
      arrange(desc(Total_CPM_Overall)) %>%
      pull(Annotation)
  }
  
  cpm_long_filtered$Annotation <- factor(cpm_long_filtered$Annotation, levels = order_annotations) # Order for vertical plot
  
  # Add granularity label to title and filename if provided
  if (!is.null(granularity_label)) {
    plot_title <- paste0(plot_title, " (", granularity_label, ")")
    output_file_name_suffix <- paste0(gsub(" ", "_", granularity_label), "_", output_file_name_suffix)
  }
  
  p <- ggplot(cpm_long_filtered, aes(x = Annotation, y = Accumulated_CPM, fill = Condition)) + # Annotation on X for vertical plot
    geom_boxplot(outlier.shape = NA) +
    geom_jitter(aes(color = Condition), width = 0.2, alpha = 0.6, size = 1) +
    scale_y_log10(
      breaks = scales::log_breaks(n = 6, base = 10),
      labels = scales::label_number(scale_cut = scales::cut_short_scale()),
      # MODIFICATION HERE: Adjust limits to extend the plot vertically
      limits = c(
        max(PSEUDOCOUNT, min(cpm_long_filtered$Accumulated_CPM[cpm_long_filtered$Accumulated_CPM > 0], na.rm = TRUE) * 0.1), # Dynamic lower bound
        max(cpm_long_filtered$Accumulated_CPM, na.rm = TRUE) * 2 # Increased upper limit for more vertical spread
      )
    ) +
    labs(
      title = plot_title,
      x = annotation_type_name,
      y = "Accumulated CPM (log10)"
    ) +
    scale_fill_brewer(palette = "Set1") +
    scale_color_brewer(palette = "Set1") +
    # NO coord_flip() for vertical boxplots
    theme_minimal(base_size = 10) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 12),
      axis.title = element_text(size = 9),
      axis.text.x = element_text(angle = 45, hjust = 1, size = 7), # Angle for x-axis labels (annotation names)
      axis.text.y = element_text(size = 8), # Y-axis is now numerical
      legend.position = "right",
      legend.title = element_text(size = 9, face = "bold"),
      panel.grid.major.x = element_blank(), # Remove vertical grid lines
      panel.grid.minor.x = element_blank(),
      axis.line = element_line(color = "black", linewidth = 0.3)
    )
  
  output_tiff <- file.path(output_dir, paste0("3_Pre_DEA_Accumulated_CPM_", gsub(" ", "_", annotation_type_name), "_", output_file_name_suffix, ".tiff"))
  output_png <- file.path(output_dir, paste0("3_Pre_DEA_Accumulated_CPM_", gsub(" ", "_", annotation_type_name), "_", output_file_name_suffix, ".png"))
  ggsave(output_tiff, plot = p, width = PLOT_WIDTH_CM, height = PLOT_HEIGHT_CM * 0.75, units = "cm", dpi = dpi, compression = "lzw")
  ggsave(output_png, plot = p, width = PLOT_WIDTH_CM, height = PLOT_HEIGHT_CM * 0.75, units = "cm", dpi = dpi)
  cat(paste0("  Accumulated CPM plot for ", annotation_type_name, 
             ifelse(!is.null(granularity_label), paste0(" (", granularity_label, ")"), ""), 
             " (", output_file_name_suffix, ") saved to: ", output_tiff, " and ", output_png, "\n"))
  return(p)
}


#' Performs enrichment analysis for a given gene list and annotation type.
#'
#' @param gene_list Character vector of significant protein accessions.
#' @param background_genes Character vector of all background protein accessions.
#' @param full_annotations_df Consolidated dataframe with all annotations.
#' @param interpro_raw_df Raw InterProScan data for method mapping.
#' @param annotation_type Type of annotation for enrichment (e.g., "GO_Terms", "EC_Numbers", "KO_Numbers", "InterPro_Methods").
#' @param comparison_name Name of the differential expression comparison.
#' @param direction Direction of change ("OverRepresented" or "UnderRepresented").
#' @param output_dir Directory to save results.
#' @param dpi Plot resolution.
#' @param analysis_prefix Prefix for output files (e.g., "Pre_DEA_FC_based" or "Post_DEA").
#' @param go_category Specific GO category to analyze (e.g., "biological_process", "molecular_function", "cellular_component"). NULL for general GO.
#' @param interpro_db Specific InterPro database to analyze (e.g., "Pfam", "PANTHER", "SUPERFAMILY", "CDD", "Gene3D", "TIGRFAMs"). NULL for general InterPro methods.
#' @return A dataframe of enrichment results, or NULL if no significant results or an error occurred.
perform_enrichment <- function(gene_list, background_genes, full_annotations_df, interpro_raw_df, annotation_type, comparison_name, direction, output_dir, dpi, analysis_prefix = "Post_DEA", go_category = NULL, interpro_db = NULL) {
  
  cat(paste0("\nPerforming enrichment for ", analysis_prefix, " ", annotation_type, 
             ifelse(!is.null(go_category), paste0(" (", go_category, ")"), ""),
             ifelse(!is.null(interpro_db), paste0(" (", interpro_db, ")"), ""),
             " - ", comparison_name, " (", direction, ")\n"))
  
  # Ensure gene_list and background_genes are unique
  gene_list <- unique(gene_list)
  background_genes <- unique(background_genes)
  
  if (length(gene_list) == 0) {
    cat("  No significant genes for this category. Skipping enrichment.\n")
    return(NULL)
  }
  
  if (length(background_genes) < 2) { # Need at least 2 for meaningful background
    cat("  Not enough background genes for enrichment analysis. Skipping.\n")
    return(NULL)
  }
  
  enrichment_results <- NULL
  tryCatch({
    if (annotation_type == "GO_Terms") {
      go_map <- full_annotations_df %>% # Use the full_annotations_df, which should have all GO rows
        filter(Protein_Accession %in% background_genes) %>%
        dplyr::select(Protein_Accession, Annotation.GO.ID, Annotation.GO.Term, Annotation.GO.Category) %>%
        filter(!is.na(Annotation.GO.ID) & Annotation.GO.ID != "") %>%
        # Ensure categories are consistently formatted and filtered
        filter(!is.na(Annotation.GO.Category) & Annotation.GO.Category %in% c("biological_process", "molecular_function", "cellular_component")) %>%
        distinct() # Keep distinct Protein_Accession-Annotation.GO.ID pairs
      
      if (!is.null(go_category)) {
        go_map <- go_map %>% filter(Annotation.GO.Category == go_category)
      }
      
      if (nrow(go_map) < 5) { # Ensure sufficient mapping data for enrichment
        cat("  Not enough valid GO annotations in the background for this type and category. Skipping.\n")
        return(NULL)
      }
      gene2term <- go_map %>% dplyr::select(Annotation.GO.ID, Protein_Accession)
      term2name <- go_map %>% dplyr::select(Annotation.GO.ID, Annotation.GO.Term) %>% distinct() %>%
        # Ensure Description is not empty for plotting by falling back to ID if term is missing
        mutate(Annotation.GO.Term = ifelse(is.na(Annotation.GO.Term) | Annotation.GO.Term == "", Annotation.GO.ID, Annotation.GO.Term))
      
      enrichment_results <- enricher(
        gene = gene_list, universe = background_genes, pAdjustMethod = "BH",
        TERM2GENE = gene2term, TERM2NAME = term2name, minGSSize = 5, maxGSSize = 500
      )
    } else if (annotation_type == "EC_Numbers") {
      ec_map <- full_annotations_df %>%
        filter(Protein_Accession %in% background_genes) %>%
        dplyr::select(Protein_Accession, Enzyme.Code, Enzyme.Name) %>%
        filter(!is.na(Enzyme.Code) & Enzyme.Code != "") %>%
        distinct()
      
      if (nrow(ec_map) < 5) { # Ensure sufficient mapping data
        cat("  Not enough valid EC number annotations in the background for this type. Skipping.\n")
        return(NULL)
      }
      gene2term <- ec_map %>% dplyr::select(Enzyme.Code, Protein_Accession)
      term2name <- ec_map %>% dplyr::select(Enzyme.Code, Enzyme.Name) %>% distinct() %>%
        mutate(Enzyme.Name = ifelse(is.na(Enzyme.Name) | Enzyme.Name == "", Enzyme.Code, Enzyme.Name))
      
      enrichment_results <- enricher(
        gene = gene_list, universe = background_genes, pAdjustMethod = "BH",
        TERM2GENE = gene2term, TERM2NAME = term2name, minGSSize = 5, maxGSSize = 500
      )
    } else if (annotation_type == "KO_Numbers") {
      ko_map <- full_annotations_df %>%
        filter(Protein_Accession %in% background_genes) %>%
        dplyr::select(Protein_Accession, KO_Number) %>%
        filter(!is.na(KO_Number) & KO_Number != "") %>%
        distinct()
      
      if (nrow(ko_map) < 5) { # Ensure sufficient mapping data
        cat("  Not enough valid KO number annotations in the background for this type. Skipping.\n")
        return(NULL)
      }
      
      # Map significant proteins to their KOs
      gene_list_kos <- ko_map %>% filter(Protein_Accession %in% gene_list) %>% pull(KO_Number) %>% unique()
      
      # Map all background proteins to their KOs
      universe_kos <- ko_map %>% pull(KO_Number) %>% unique()
      
      if (length(gene_list_kos) == 0 || length(universe_kos) < 5) {
        cat("  Not enough KOs for enrichment after filtering. Skipping KEGG enrichment.\n")
        return(NULL)
      }
      
      enrichment_results <- enrichKEGG(
        gene = gene_list_kos, universe = universe_kos, organism = 'ko', pAdjustMethod = "BH",
        minGSSize = 5, maxGSSize = 500
      )
    } else if (annotation_type == "InterPro_Methods") {
      # For InterPro enrichment, use interpro_raw_df to get individual Analysis_Method entries
      interpro_method_map <- interpro_raw_df %>%
        filter(Protein_Accession %in% background_genes & Status == "T") %>%
        dplyr::select(Protein_Accession, Analysis_Method, Signature_Description) %>% # Include Signature_Description for better TERM2NAME
        filter(!is.na(Analysis_Method) & Analysis_Method != "") %>% # Analysis_Method should not be NA or empty
        distinct() # Keep distinct Protein_Accession-Analysis_Method pairs
      
      if (!is.null(interpro_db)) {
        interpro_method_map <- interpro_method_map %>% filter(Analysis_Method == interpro_db)
      }
      
      if (nrow(interpro_method_map) < 5) { # Ensure sufficient mapping data
        cat("  Not enough valid InterPro Method annotations in the background for this type and database. Skipping.\n")
        return(NULL)
      }
      
      gene2term <- interpro_method_map %>% dplyr::select(Analysis_Method, Protein_Accession)
      
      # Create TERM2NAME using Analysis_Method as ID and Signature_Description as Term where available, otherwise use ID
      term2name <- interpro_method_map %>% 
        distinct(Analysis_Method, Signature_Description) %>%
        dplyr::rename(ID = Analysis_Method) %>% 
        # For term2name, if Signature_Description is available and not NA/empty, use it. Otherwise, use ID.
        mutate(Term = ifelse(!is.na(Signature_Description) & Signature_Description != "", Signature_Description, ID)) %>%
        dplyr::select(ID, Term) %>%
        distinct() # Ensure unique ID-Term pairs
      
      # Fallback: if after all filtering, term2name is empty or has issues
      if (nrow(term2name) == 0) {
        # Fallback to just Analysis_Method as term name if descriptions are all missing
        term2name <- interpro_method_map %>% distinct(Analysis_Method) %>%
          dplyr::rename(ID = Analysis_Method) %>% mutate(Term = ID)
      }
      
      if (ncol(gene2term) != 2 || ncol(term2name) != 2 || nrow(gene2term) == 0 || nrow(term2name) == 0 ) {
        cat("  Malformed gene2term or term2name for InterPro_Methods (dimensions/empty). Skipping enrichment.\n")
        return(NULL)
      }
      
      enrichment_results <- enricher(
        gene = gene_list, universe = background_genes, pAdjustMethod = "BH",
        TERM2GENE = gene2term, TERM2NAME = term2name, minGSSize = 5, maxGSSize = 500
      )
    } else {
      cat("  Enrichment for '", annotation_type, "' is not yet implemented or recognized. Skipping.\n")
      return(NULL)
    }
    
    if (!is.null(enrichment_results) && nrow(enrichment_results@result) > 0) {
      cat(paste0("  Found ", nrow(enrichment_results@result), " significant enrichments for ", annotation_type, ".\n"))
      
      enrich_df <- enrichment_results@result
      
      # Ensure Description is not empty for plotting by falling back to ID if term is missing
      if ("Description" %in% colnames(enrich_df)) {
        enrich_df$Description <- ifelse(is.na(enrich_df$Description) | enrich_df$Description == "", enrich_df$ID, enrich_df$Description)
      } else {
        enrich_df$Description <- enrich_df$ID # Fallback if Description column is missing
      }
      
      enrich_df <- enrich_df %>%
        dplyr::rename(Annotation.ID = ID) %>%
        mutate(
          Annotation.Type = annotation_type,
          Comparison = comparison_name,
          Direction = direction,
          GO.Category = ifelse(!is.null(go_category), go_category, NA_character_), # Add GO category
          InterPro.DB = ifelse(!is.null(interpro_db), interpro_db, NA_character_) # Add InterPro DB
        ) %>%
        dplyr::select(Annotation.Type, Comparison, Direction, GO.Category, InterPro.DB, Annotation.ID, Description, everything())
      
      if (nrow(enrich_df) == 0) {
        cat("  Warning: enrich_df became empty after processing (e.g., due to filtering out all descriptions). Not writing CSV or plotting.\n")
        return(NULL)
      }
      
      file_suffix_enrich <- ""
      if (!is.null(go_category)) {
        file_suffix_enrich <- paste0("_", gsub(" ", "_", go_category))
      } else if (!is.null(interpro_db)) {
        file_suffix_enrich <- paste0("_", gsub(" ", "_", interpro_db))
      }
      
      output_enrich_csv <- file.path(output_dir, paste0(analysis_prefix, "_Enrichment_", annotation_type, file_suffix_enrich, "_", comparison_name, "_", direction, ".csv"))
      write.csv(enrich_df, output_enrich_csv, row.names = FALSE)
      cat("  Enrichment results saved to:", output_enrich_csv, "\n")
      
      # Prepare for plotting, ensuring 'Description' is available
      if ("Description" %in% colnames(enrichment_results@result)) {
        enrichment_results@result$Description <- ifelse(is.na(enrichment_results@result$Description) | enrichment_results@result$Description == "", enrichment_results@result$ID, enrichment_results@result$Description)
      } else {
        enrichment_results@result$Description <- enrichment_results@result$ID
      }
      
      tryCatch({
        num_categories_to_show <- min(15, nrow(enrichment_results@result))
        if (num_categories_to_show == 0) {
          cat("  Warning: No categories left to plot after filtering or processing for dotplot. Skipping.\n")
          return(NULL)
        }
        
        plot_data <- as.data.frame(enrichment_results@result) %>%
          arrange(p.adjust, desc(Count)) %>%
          head(num_categories_to_show) %>%
          mutate(Description = forcats::fct_reorder(Description, Count))
        
        p_enrichment_dotplot <- ggplot(plot_data, aes(x = Count, y = Description, size = Count, color = p.adjust)) +
          geom_point() +
          scale_color_gradient(low = "red", high = "blue", name = "Adjusted P-value") +
          labs(
            title = paste0("Enriched ", annotation_type, 
                           ifelse(!is.null(go_category), paste0(" (", go_category, ")"), ""),
                           ifelse(!is.null(interpro_db), paste0(" (", interpro_db, ")"), ""),
                           " (", comparison_name, ": ", direction, ")"),
            x = "Gene Count",
            y = ""
          ) +
          theme_minimal(base_size = 10) +
          theme(
            plot.title = element_text(hjust = 0.5, face = "bold", size = 12),
            axis.title = element_text(size = 9),
            axis.text.y = element_text(size = 8),
            axis.text.x = element_text(size = 8),
            legend.position = "right"
          )
        
        output_enrich_tiff <- file.path(output_dir, paste0(analysis_prefix, "_Enrichment_DotPlot_", annotation_type, file_suffix_enrich, "_", comparison_name, "_", direction, ".tiff"))
        output_enrich_png <- file.path(output_dir, paste0(analysis_prefix, "_Enrichment_DotPlot_", annotation_type, file_suffix_enrich, "_", comparison_name, "_", direction, ".png"))
        
        ggsave(output_enrich_tiff, plot = p_enrichment_dotplot, width = PLOT_WIDTH_CM * 0.7, height = PLOT_HEIGHT_CM * 0.7, units = "cm", dpi = dpi, compression = "lzw")
        ggsave(output_png, plot = p_enrichment_dotplot, width = PLOT_WIDTH_CM * 0.7, height = PLOT_HEIGHT_CM * 0.7, units = "cm", dpi = dpi)
        cat("  Enrichment DotPlot saved to: ", output_enrich_tiff, " and ", output_enrich_png, "\n")
      }, error = function(e) {
        cat("  Warning: Could not generate dotplot for", annotation_type, 
            ifelse(!is.null(go_category), paste0(" (", go_category, ")"), ""),
            ifelse(!is.null(interpro_db), paste0(" (", interpro_db, ")"), ""),
            " - ", comparison_name, " (", direction, "). Error:", e$message, "\n")
      })
      return(enrich_df)
    } else {
      cat("  No significant enrichments found for this category.\n")
    }
  }, error = function(e) {
    cat("  Error during enrichment for", annotation_type, 
        ifelse(!is.null(go_category), paste0(" (", go_category, ")"), ""),
        ifelse(!is.null(interpro_db), paste0(" (", interpro_db, ")"), ""),
        ":", e$message, "\n")
    return(NULL)
  })
}

# ==============================================================================
# 3. Main Analysis Workflow 🚀
# ==============================================================================

main_analysis_workflow <- function() {
  cat("\n==============================================================================")
  cat("\nStarting comprehensive annotation and differential expression analysis...")
  cat("\n==============================================================================\n")
  
  # Create output directory if it doesn't exist
  if (!dir.exists(OUTPUT_RESULTS_DIR)) {
    dir.create(OUTPUT_RESULTS_DIR, recursive = TRUE)
    cat("Output folder created:", OUTPUT_RESULTS_DIR, "\n")
  } else {
    cat("Using existing output folder:", OUTPUT_RESULTS_DIR, "\n")
  }
  
  # --- Step 1: Load and Consolidate All Annotation Data ---
  cat("\n--- Consolidating all annotation data ---\n")
  
  # Load main EC abundance file (source of Protein_Accession and CPMs)
  cat("Loading main data from:", FILE_PATHS$ec_abundance, "...\n")
  if (!file.exists(FILE_PATHS$ec_abundance)) {
    stop(paste("Error: Main EC abundance file DOES NOT exist at path:", FILE_PATHS$ec_abundance))
  }
  main_data_df <- readr::read_delim(FILE_PATHS$ec_abundance, delim = DELIMITERS$ec_abundance, col_names = TRUE, show_col_types = FALSE, col_types = cols(.default = col_character()))
  
  if (PROTEIN_ID_COL_NAME_MAIN %in% colnames(main_data_df)) {
    main_data_df <- main_data_df %>% dplyr::rename(Protein_Accession = !!sym(PROTEIN_ID_COL_NAME_MAIN))
    cat(paste0("Column '", PROTEIN_ID_COL_NAME_MAIN, "' renamed to 'Protein_Accession'.\n"))
  } else {
    stop(paste0("Error: Column '", PROTEIN_ID_COL_NAME_MAIN, "' (for Protein_Accession) DOES NOT exist in ", basename(FILE_PATHS$ec_abundance), ". Columns found: ", paste(colnames(main_data_df), collapse = ", ")))
  }
  
  if ("MAG" %in% colnames(main_data_df)) {
    main_data_df <- main_data_df %>% dplyr::select(-any_of("MAG"))
    cat("Column 'MAG' removed.\n")
  } else {
    warning("Warning: Column 'MAG' was not found to remove.")
  }
  
  cpm_cols_names_actual <- colnames(main_data_df)[str_detect(colnames(main_data_df), CPM_COLS_PATTERN)]
  if (length(cpm_cols_names_actual) == 0) {
    stop("Error: No CPM columns found in the main file. Check your column names and pattern.")
  } else {
    cat(paste0("Identified CPM columns: ", paste(cpm_cols_names_actual, collapse = ", "), "\n"))
  }
  
  main_data_df <- main_data_df %>%
    mutate(across(all_of(cpm_cols_names_actual), ~na_if(., "-"))) %>%
    mutate(across(all_of(cpm_cols_names_actual), as.numeric))
  
  if (any(!sapply(main_data_df %>% dplyr::select(all_of(cpm_cols_names_actual)), is.numeric))) {
    stop("Error: Some CPM columns could not be converted to numeric after cleaning. Please check your data.")
  } else {
    cat("CPM columns successfully converted to numeric in main_data_df.\n")
  }
  
  if ("EC_number" %in% colnames(main_data_df) && "product" %in% colnames(main_data_df)) {
    main_data_df <- main_data_df %>%
      dplyr::rename(Enzyme.Code = EC_number, Enzyme.Name = product) %>%
      mutate(across(where(is.character), ~na_if(., ""))) %>%
      distinct(Protein_Accession, .keep_all = TRUE)
  } else {
    stop("Error: 'EC_number' or 'product' columns not found in the main file. Check their exact names.")
  }
  cat("Main data loaded and preprocessed.\n")
  
  # Load OmicBox_annotation.csv (for GO and descriptions)
  cat("Loading OmicBox annotations from:", FILE_PATHS$omicbox_annotation, "...\n")
  if (!file.exists(FILE_PATHS$omicbox_annotation)) {
    stop(paste("Error: OmicBox annotation file DOES NOT exist at path:", FILE_PATHS$omicbox_annotation))
  }
  omicbox_raw_df <- readr::read_delim(FILE_PATHS$omicbox_annotation, delim = DELIMITERS$omicbox_annotation, col_names = TRUE, show_col_types = FALSE, col_types = cols(.default = col_character()))
  
  # *** MODIFIED GO PROCESSING ***
  omicbox_processed_df <- omicbox_raw_df %>%
    dplyr::rename(Protein_Accession = Sequence_Name, Sequence.Description = Sequence_Description,
                  `Blast.Top.Hit.Taxonomy.Name` = `Blast_Top_Hit_Taxonomy Name`,
                  Annotation.GO.ID = Annotation_GO_ID, Annotation.GO.Term = Annotation_GO_Term,
                  Annotation.GO.Category = Annotation_GO_Category) %>%
    dplyr::select(Protein_Accession, any_of(c("Sequence.Description", "Blast.Top.Hit.Taxonomy.Name", "Annotation.GO.ID", "Annotation.GO.Term", "Annotation.GO.Category"))) %>%
    mutate(
      Annotation.GO.ID = str_extract(Annotation.GO.ID, "(go|GO):\\d{7}"),
      Annotation.GO.ID = str_to_upper(Annotation.GO.ID),
      # Robust standardization for Annotation.GO.Category
      Annotation.GO.Category = tolower(trimws(Annotation.GO.Category)), # Trim spaces first
      Annotation.GO.Category = case_when(
        grepl("biological process", Annotation.GO.Category) ~ "biological_process",
        grepl("molecular function", Annotation.GO.Category) ~ "molecular_function",
        grepl("cellular component", Annotation.GO.Category) ~ "cellular_component",
        TRUE ~ NA_character_ # Assign NA if it doesn't match expected categories
      )
    ) %>%
    mutate(across(where(is.character), ~na_if(., "")))
  # We need to keep all rows for correct GO summarization and filtering.
  
  cat("OmicBox annotation data cleaned and standardized.\n")
  
  # Load user_ko.txt (PRIMARY KO SOURCE)
  cat("Loading KO annotations from:", FILE_PATHS$ko_annotations, "...\n")
  if (!file.exists(FILE_PATHS$ko_annotations)) {
    stop(paste("Error: KO annotation file DOES NOT exist at path:", FILE_PATHS$ko_annotations))
  }
  ko_primary_source_df <- readr::read_lines(FILE_PATHS$ko_annotations) %>%
    tibble(raw_line = .) %>%
    filter(!startsWith(trimws(raw_line), "#")) %>%
    filter(nchar(trimws(raw_line)) > 0) %>%
    separate_wider_delim(raw_line, delim = DELIMITERS$ko_annotations, names = c("Protein_Accession", "KO_Number"), too_few = "align_start", too_many = "merge") %>%
    mutate(
      Protein_Accession = trimws(Protein_Accession),
      KO_Number = trimws(KO_Number)
    ) %>%
    filter(!is.na(Protein_Accession) & Protein_Accession != "") %>%
    mutate(across(where(is.character), ~na_if(., ""))) %>%
    distinct(Protein_Accession, .keep_all = TRUE) %>%
    mutate(KO_Name = NA_character_)
  cat("Primary KO source file loaded and processed.\n")
  
  # Load raw InterPro data for specific InterPro annotations
  cat("Loading raw InterProScan data from:", FILE_PATHS$interproscan_full, "...\n")
  if (!file.exists(FILE_PATHS$interproscan_full)) {
    stop(paste("Error: InterProScan file DOES NOT exist at path:", FILE_PATHS$interproscan_full))
  }
  interpro_raw_df <- readr::read_delim(
    FILE_PATHS$interproscan_full,
    delim = DELIMITERS$interproscan_full,
    col_names = c("Protein_Accession", "MD5_Digest", "Sequence_Length", "Analysis_Method", "Signature_Accession", "Signature_Description", "Start_Location", "End_Location", "E_value", "Status", "Date", "InterPro_Accession", "InterPro_Description", "GO_Annotations", "Pathways"),
    col_types = cols(.default = col_character())
  )
  
  # ADDED: Robustly convert "-" or "" to NA for InterPro annotation columns
  interpro_raw_df <- interpro_raw_df %>%
    mutate(
      Signature_Accession = na_if(Signature_Accession, "-"),
      Signature_Accession = na_if(Signature_Accession, ""),
      Signature_Description = na_if(Signature_Description, "-"),
      Signature_Description = na_if(Signature_Description, ""),
      InterPro_Description = na_if(InterPro_Description, "-"),
      InterPro_Description = na_if(InterPro_Description, "")
    )
  cat("InterProScan raw data loaded and '-, \"\" ' values converted to NA.\n")
  
  # This interpro_summary_df is for a consolidated view, not for granular plotting
  interpro_summary_df <- interpro_raw_df %>%
    filter(Status == "T") %>%
    group_by(Protein_Accession) %>%
    summarise(
      InterPro_Domains = paste(unique(na.omit(InterPro_Description)), collapse = " | "),
      Pfam_Domains = paste(unique(na.omit(Signature_Description[Analysis_Method == "Pfam"])), collapse = " | "),
      PANTHER_Domains = paste(unique(na.omit(Signature_Description[Analysis_Method == "PANTHER"])), collapse = " | "),
      CDD_Domains = paste(unique(na.omit(Signature_Description[Analysis_Method == "CDD"])), collapse = " | "),
      NCBIfam_Domains = paste(unique(na.omit(Signature_Description[Analysis_Method == "NCBIfam"])), collapse = " | "),
      ProSiteProfiles_Domains = paste(unique(na.omit(Signature_Description[Analysis_Method == "ProSiteProfiles"])), collapse = " | "),
      SMART_Domains = paste(unique(na.omit(Signature_Description[Analysis_Method == "SMART"])), collapse = " | "),
      SUPERFAMILY_Domains = paste(unique(na.omit(Signature_Description[Analysis_Method == "SUPERFAMILY"])), collapse = " | "),
      Gene3D_Domains = paste(unique(na.omit(Signature_Description[Analysis_Method == "Gene3D"])), collapse = " | "),
      TIGRFAMs_Domains = paste(unique(na.omit(Signature_Description[Analysis_Method == "TIGRFAMs"])), collapse = " | "),
      N_Annotations_InterPro = n(),
      .groups = 'drop'
    ) %>%
    mutate(across(where(is.character), ~na_if(., ""))) %>%
    mutate(InterPro_Domains = ifelse(InterPro_Domains == "-", NA_character_, InterPro_Domains))
  cat("InterProScan summary data prepared.\n")
  
  # Build the FINAL full_annotations_df
  # For GO, we need to join back omicbox_processed_df *without* the distinct(Protein_Accession) step,
  # or ensure it's handled correctly for many-to-one protein-GO relationships.
  # Let's create `full_annotations_with_all_gos` which allows multiple GOs per protein.
  
  full_annotations_with_all_gos <- main_data_df %>%
    left_join(omicbox_processed_df, by = "Protein_Accession", relationship = "many-to-many") %>% # Use many-to-many for GO
    left_join(ko_primary_source_df, by = "Protein_Accession") %>%
    # Note: interpro_summary_df has one row per protein. For granular InterPro plots, we'll use interpro_raw_df
    left_join(interpro_summary_df, by = "Protein_Accession") %>% 
    dplyr::select(
      Protein_Accession, all_of(cpm_cols_names_actual),
      Enzyme.Code, Enzyme.Name, Annotation.GO.ID, Annotation.GO.Term, Annotation.GO.Category,
      KO_Number, KO_Name, Sequence.Description, Blast.Top.Hit.Taxonomy.Name,
      InterPro_Domains, Pfam_Domains, PANTHER_Domains, CDD_Domains, NCBIfam_Domains,
      ProSiteProfiles_Domains, SMART_Domains, SUPERFAMILY_Domains, Gene3D_Domains, TIGRFAMs_Domains,
      N_Annotations_InterPro
    ) %>%
    mutate(across(where(is.character), ~na_if(., "")))
  
  # For limma and general background, we still need a version where each protein is unique.
  # This `full_annotations_df_for_limma_and_background` will be used for overall DE analysis and background gene lists.
  full_annotations_df_for_limma_and_background <- full_annotations_with_all_gos %>%
    distinct(Protein_Accession, .keep_all = TRUE) # This version is safe for limma
  
  output_full_annotations_csv <- file.path(OUTPUT_RESULTS_DIR, "Full_Consolidated_Annotations_Tenebrio_Molitor.csv")
  write.csv(full_annotations_with_all_gos, output_full_annotations_csv, row.names = FALSE) # Save the one with all GO rows
  cat("Consolidated annotation file saved to:", output_full_annotations_csv, "\n")
  cat("Please review 'Full_Consolidated_Annotations_Tenebrio_Molitor.csv' to confirm GO, EC, and KO categories are correct.\n")
  
  # --- Step 2: Run Annotation Abundance and Cumulative Counts Analysis ---
  cat("\n--- Running InterProScan annotation abundance and cumulative counts analysis ---\n")
  
  # Plot Abundance of InterProScan Methods (overall)
  interpro_with_cpm_for_plot <- interpro_raw_df %>%
    filter(Status == "T") %>% # Only include 'True' matches
    dplyr::select(Protein_Accession, Analysis_Method) %>%
    distinct() %>% # Ensure unique protein-method pairs
    left_join(
      main_data_df %>% # Join to main_data_df for CPMs
        dplyr::select(Protein_Accession, starts_with(c("C_", "LD_", "PET_"))),
      by = "Protein_Accession"
    ) %>%
    filter(!is.na(Analysis_Method))
  
  p_abundance_method <- NULL
  if (nrow(interpro_with_cpm_for_plot) > 0) {
    interpro_method_abundance_df <- interpro_with_cpm_for_plot %>%
      group_by(Analysis_Method) %>%
      summarise(across(starts_with(c("C_", "LD_", "PET_")), \(x) sum(x, na.rm = TRUE)), .groups = 'drop') %>% # Fix dplyr across() warning
      dplyr::rename(Method = Analysis_Method)
    
    p_abundance_method <- plot_abundance_by_method(
      data = interpro_method_abundance_df, output_dir = OUTPUT_RESULTS_DIR, file_suffix = "Overall_InterProMethods",
      width_cm = PLOT_WIDTH_CM / 2, height_cm = PLOT_HEIGHT_CM / 2, dpi = DPI
    )
  } else {
    cat("No InterProScan data with CPMs available for plotting overall abundance.\n")
  }
  
  # Plot Cumulative Unique Annotations for main types
  
  # 1) Overall EC Numbers
  p_cumulative_ec_overall <- plot_cumulative_annotations(
    data = full_annotations_df_for_limma_and_background,
    annotation_col = "Enzyme.Code",
    annotation_type_name = "EC Numbers",
    output_dir = OUTPUT_RESULTS_DIR, file_suffix = "EC_Numbers", # Simplified suffix
    width_cm = PLOT_WIDTH_CM / 2, height_cm = PLOT_HEIGHT_CM / 2, dpi = DPI
  )
  
  # 2) Overall KO Numbers
  p_cumulative_ko_overall <- plot_cumulative_annotations(
    data = full_annotations_df_for_limma_and_background,
    annotation_col = "KO_Number",
    annotation_type_name = "KO Numbers",
    output_dir = OUTPUT_RESULTS_DIR, file_suffix = "KO_Numbers", # Simplified suffix
    width_cm = PLOT_WIDTH_CM / 2, height_cm = PLOT_HEIGHT_CM / 2, dpi = DPI
  )
  
  # 3) GO Terms (Biological Process, Molecular Function, Cellular Component)
  go_categories_for_plots <- c("biological_process", "molecular_function", "cellular_component")
  for (cat_name in go_categories_for_plots) {
    cat(paste0("Generating cumulative unique annotations plot for GO: ", str_to_title(gsub("_", " ", cat_name)), "...\n"))
    plot_cumulative_annotations(
      data = full_annotations_with_all_gos %>% filter(Annotation.GO.Category == cat_name),
      annotation_col = "Annotation.GO.ID",
      annotation_type_name = "GO Terms",
      output_dir = OUTPUT_RESULTS_DIR, file_suffix = paste0("GO_", str_to_title(gsub("_", "", cat_name))), # Simplified suffix
      width_cm = PLOT_WIDTH_CM / 2, height_cm = PLOT_HEIGHT_CM / 2, dpi = DPI
    )
  }
  
  # 4) InterPro Methods (Specific Databases)
  interpro_databases_for_plots <- c("Pfam", "PANTHER", "SUPERFAMILY", "CDD", "Gene3D", "TIGRFAMs")
  for (db_name in interpro_databases_for_plots) {
    cat(paste0("Generating cumulative unique annotations plot for InterPro: ", db_name, "...\n"))
    plot_cumulative_annotations(
      data = interpro_raw_df %>% filter(Status == "T", Analysis_Method == db_name),
      annotation_col = "Signature_Accession", # Use Signature_Accession for unique domain types
      annotation_type_name = "InterPro Methods",
      output_dir = OUTPUT_RESULTS_DIR, file_suffix = paste0("InterPro_", db_name), # Simplified suffix
      width_cm = PLOT_WIDTH_CM / 2, height_cm = PLOT_HEIGHT_CM / 2, dpi = DPI
    )
  }
  
  # REMOVED: The section that was causing the 'p_cumulative_interpro_overall' error.
  # This section was a remnant from a previous logic and is no longer needed.
  
  cat("\nAnnotation abundance and cumulative counts analysis completed.\n")
  
  # --- Step 3: Pre-DEA Enrichment: Accumulated CPM Boxplots ---
  cat("\n--- Performing Pre-DEA Enrichment: Accumulated CPM Boxplots per Annotation Type and Condition (Granular) ---\n")
  
  conditions_for_top15_plots <- c("C", "LD", "PET")
  
  # GO Terms - granular boxplots
  go_categories_for_plots <- c("biological_process", "molecular_function", "cellular_component")
  for (cat_name in go_categories_for_plots) {
    # Use full_annotations_with_all_gos to ensure all GO rows are considered
    df_for_go_plot <- full_annotations_with_all_gos %>%
      filter(Annotation.GO.Category == cat_name)
    
    # Global Top 15 for this GO category
    plot_accumulated_cpm(df_for_go_plot, cpm_cols_names_actual, "Annotation.GO.ID", "GO Terms", 
                         OUTPUT_RESULTS_DIR, DPI, plot_type = "global_top15", granularity_label = paste0("GO: ", str_to_title(gsub("_", " ", cat_name))))
    
    # Top 15 per condition for this GO category
    for (cond in conditions_for_top15_plots) {
      plot_accumulated_cpm(df_for_go_plot, cpm_cols_names_actual, "Annotation.GO.ID", "GO Terms", 
                           OUTPUT_RESULTS_DIR, DPI, plot_type = cond, granularity_label = paste0("GO: ", str_to_title(gsub("_", " ", cat_name))))
    }
  }
  
  # InterPro Methods - granular boxplots
  interpro_databases_for_plots <- c("Pfam", "PANTHER", "SUPERFAMILY", "CDD", "Gene3D", "TIGRFAMs")
  for (db_name in interpro_databases_for_plots) {
    # *** MODIFIED INTERPRO PLOT DATA PREPARATION ***
    df_for_interpro_plot <- interpro_raw_df %>%
      filter(Status == "T", Analysis_Method == db_name) %>% # Filter for specific DB and valid status
      dplyr::select(Protein_Accession, Signature_Description) %>% # Select Protein and the specific description
      distinct() %>% # Ensure unique protein-signature description pairs
      left_join(
        main_data_df %>% # Join CPMs from main_data_df
          dplyr::select(Protein_Accession, all_of(cpm_cols_names_actual)),
        by = "Protein_Accession"
      ) %>%
      filter(!is.na(Signature_Description)) # Only keep rows with a valid signature description (NA handled upstream)
    
    # Global Top 15 for this InterPro database
    plot_accumulated_cpm(df_for_interpro_plot, cpm_cols_names_actual, "Signature_Description", "InterPro Methods", # Plot by Signature_Description
                         OUTPUT_RESULTS_DIR, DPI, plot_type = "global_top15", granularity_label = db_name)
    
    # Top 15 per condition for this InterPro database
    for (cond in conditions_for_top15_plots) {
      plot_accumulated_cpm(df_for_interpro_plot, cpm_cols_names_actual, "Signature_Description", "InterPro Methods", # Plot by Signature_Description
                           OUTPUT_RESULTS_DIR, DPI, plot_type = cond, granularity_label = db_name)
    }
  }
  
  # EC Numbers - overall and per condition plots
  # Use full_annotations_df_for_limma_and_background for EC/KO to ensure one protein per row, as that's how EC/KO are generally structured.
  plot_accumulated_cpm(full_annotations_df_for_limma_and_background, cpm_cols_names_actual, "Enzyme.Code", "EC Numbers", OUTPUT_RESULTS_DIR, DPI, plot_type = "global_top15")
  for (cond in conditions_for_top15_plots) {
    plot_accumulated_cpm(full_annotations_df_for_limma_and_background, cpm_cols_names_actual, "Enzyme.Code", "EC Numbers", OUTPUT_RESULTS_DIR, DPI, plot_type = cond)
  }
  
  # KO Numbers - overall and per condition plots
  plot_accumulated_cpm(full_annotations_df_for_limma_and_background, cpm_cols_names_actual, "KO_Number", "KO Numbers", OUTPUT_RESULTS_DIR, DPI, plot_type = "global_top15")
  for (cond in conditions_for_top15_plots) {
    plot_accumulated_cpm(full_annotations_df_for_limma_and_background, cpm_cols_names_actual, "KO_Number", "KO Numbers", OUTPUT_RESULTS_DIR, DPI, plot_type = cond)
  }
  
  cat("\n--- Pre-DEA Accumulated CPMs analysis completed. ---\n")
  
  # --- Step 4: Pre-DEA Fold Change Calculation and Enrichment (NEW STEP) ---
  cat("\n--- Performing Pre-DEA Fold Change Calculation and Enrichment ---\n")
  
  # Calculate mean CPMs per protein per condition
  mean_cpm_df <- full_annotations_df_for_limma_and_background %>% # Use distinct version for mean CPM calculation
    dplyr::select(Protein_Accession, all_of(cpm_cols_names_actual)) %>%
    pivot_longer(
      cols = all_of(cpm_cols_names_actual),
      names_to = "Sample",
      values_to = "CPM_Value"
    ) %>%
    mutate(Condition = sub("_[0-9]+K", "", Sample)) %>%
    group_by(Protein_Accession, Condition) %>%
    summarise(Mean_CPM = mean(CPM_Value, na.rm = TRUE), .groups = 'drop') %>%
    pivot_wider(names_from = Condition, values_from = Mean_CPM, values_fill = 0)
  
  # Comparisons for Pre-DEA FC
  pre_dea_comparisons <- list(
    c("C", "LD"),
    c("C", "PET"),
    c("LD", "PET")
  )
  
  pre_dea_significant_proteins <- list()
  background_proteins <- unique(full_annotations_df_for_limma_and_background$Protein_Accession) # Same background as Post-DEA
  
  for (comp_pair in pre_dea_comparisons) {
    cond1 <- comp_pair[1]
    cond2 <- comp_pair[2]
    comp_name <- paste0(cond1, "_vs_", cond2)
    cat(paste0("  - Calculating Pre-DEA FC for ", comp_name, "...\n"))
    
    fc_df <- mean_cpm_df %>%
      mutate(
        # Add a pseudocount for log2 calculation if any mean CPM is 0
        Mean_CPM_cond1 = .data[[cond1]] + PSEUDOCOUNT, # Small pseudocount
        Mean_CPM_cond2 = .data[[cond2]] + PSEUDOCOUNT,
        Log2FC = log2(Mean_CPM_cond1 / Mean_CPM_cond2)
      ) %>%
      dplyr::select(Protein_Accession, Log2FC)
    
    # Categorize based on pre-defined threshold
    over_represented_proteins <- fc_df %>%
      filter(Log2FC > PRE_DEA_LOG2FC_THRESHOLD) %>%
      pull(Protein_Accession)
    
    under_represented_proteins <- fc_df %>%
      filter(Log2FC < -PRE_DEA_LOG2FC_THRESHOLD) %>%
      pull(Protein_Accession)
    
    pre_dea_significant_proteins[[comp_name]]$over <- over_represented_proteins
    pre_dea_significant_proteins[[comp_name]]$under <- under_represented_proteins
    
    cat(paste0("  Pre-DEA FC based ", comp_name, ":\n    Over-represented proteins (log2FC > ", PRE_DEA_LOG2FC_THRESHOLD, "): ", length(over_represented_proteins), "\n    Under-represented proteins (log2FC < -", PRE_DEA_LOG2FC_THRESHOLD, "): ", length(under_represented_proteins), "\n"))
    
    # Perform enrichment for Pre-DEA FC based lists
    for (direction in c("OverRepresented", "UnderRepresented")) {
      gene_list_pre_dea <- if (direction == "OverRepresented") over_represented_proteins else under_represented_proteins
      
      # For GO Terms, iterate through categories
      go_categories <- c("biological_process", "molecular_function", "cellular_component")
      for (cat in go_categories) {
        perform_enrichment(
          gene_list = gene_list_pre_dea,
          background_genes = background_proteins,
          full_annotations_df = full_annotations_with_all_gos, # Use all GO rows for enrichment mapping
          interpro_raw_df = interpro_raw_df,
          annotation_type = "GO_Terms",
          comparison_name = comp_name,
          direction = direction,
          output_dir = OUTPUT_RESULTS_DIR,
          dpi = DPI,
          analysis_prefix = "Pre_DEA_FC_based",
          go_category = cat # Pass the specific GO category
        )
      }
      
      # For InterPro Methods, iterate through specific databases
      interpro_databases <- c("Pfam", "PANTHER", "SUPERFAMILY", "CDD", "Gene3D", "TIGRFAMs")
      for (db in interpro_databases) {
        perform_enrichment(
          gene_list = gene_list_pre_dea,
          background_genes = background_proteins,
          full_annotations_df = full_annotations_df_for_limma_and_background, # Still use distinct for this, as interpro_raw_df is passed too
          interpro_raw_df = interpro_raw_df, # Raw interpro data is crucial for granular enrichment
          annotation_type = "InterPro_Methods",
          comparison_name = comp_name,
          direction = direction,
          output_dir = OUTPUT_RESULTS_DIR,
          dpi = DPI,
          analysis_prefix = "Pre_DEA_FC_based",
          interpro_db = db # Pass the specific InterPro database
        )
      }
      
      # For EC Numbers and KO Numbers (general analysis)
      for (enrich_type in c("EC_Numbers", "KO_Numbers")) {
        perform_enrichment(
          gene_list = gene_list_pre_dea,
          background_genes = background_proteins,
          full_annotations_df = full_annotations_df_for_limma_and_background,
          interpro_raw_df = interpro_raw_df, # Keep, as it's a general parameter
          annotation_type = enrich_type,
          comparison_name = comp_name,
          direction = direction,
          output_dir = OUTPUT_RESULTS_DIR,
          dpi = DPI,
          analysis_prefix = "Pre_DEA_FC_based"
        )
      }
    }
  }
  cat("\n--- Pre-DEA Fold Change and Enrichment analysis completed. ---\n")
  
  
  # --- Step 5: Prepare data for limma (Directly on log2 CPM) ---
  cat("\n--- Preparing data for differential expression analysis with limma (log2 CPM) ---\n")
  # Use the distinct version for limma, as it expects one row per gene/protein
  cpm_cols_matrix <- as.matrix(full_annotations_df_for_limma_and_background %>% dplyr::select(all_of(cpm_cols_names_actual)))
  cpm_cols_matrix[is.na(cpm_cols_matrix)] <- 0
  
  keep_expressed_anywhere <- rowSums(cpm_cols_matrix > 0) > 0
  cpm_cols_final_filtered <- cpm_cols_matrix[keep_expressed_anywhere, ]
  
  if (nrow(cpm_cols_final_filtered) == 0) {
    stop("After removing all-zero rows, no proteins remain. Check your data, they might all be zero.")
  }
  
  protein_accessions_filtered <- full_annotations_df_for_limma_and_background$Protein_Accession[keep_expressed_anywhere]
  rownames(cpm_cols_final_filtered) <- protein_accessions_filtered
  
  logCPM_values <- limma::voom(cpm_cols_final_filtered, plot=FALSE, lib.size=colSums(cpm_cols_final_filtered))$E # Using voom for better variance stabilization
  
  conditions <- factor(c(rep("C", 4), rep("LD", 4), rep("PET", 4)), levels = c("C", "LD", "PET"))
  design <- model.matrix(~0 + conditions)
  colnames(design) <- levels(conditions)
  
  contr.matrix <- makeContrasts(C_vs_LD = C - LD, C_vs_PET = C - PET, LD_vs_PET = LD - PET, levels = design)
  
  cat("\nPerforming differential expression analysis with limma on log2 transformed CPMs...\n")
  fit <- lmFit(logCPM_values, design)
  fit2 <- contrasts.fit(fit, contr.matrix)
  fit2 <- eBayes(fit2)
  
  comparisons <- colnames(contr.matrix)
  all_de_results <- list()
  for (comp in comparisons) {
    cat("  - Analyzing comparison:", comp, "\n")
    top_table <- topTable(fit2, coef = comp, number = Inf, adjust.method = "BH") %>%
      rownames_to_column("Protein_Accession") %>%
      mutate(Comparison = comp)
    
    # Join with the distinct version of full_annotations_df
    top_table_annotated <- top_table %>%
      left_join(dplyr::select(full_annotations_df_for_limma_and_background, Protein_Accession, Enzyme.Code, Enzyme.Name,
                              InterPro_Domains, Pfam_Domains, Annotation.GO.ID, Annotation.GO.Term,
                              Sequence.Description, Annotation.GO.Category, KO_Number, KO_Name,
                              SUPERFAMILY_Domains, Gene3D_Domains, TIGRFAMs_Domains),
                by = "Protein_Accession")
    all_de_results[[comp]] <- top_table_annotated
    output_de_csv <- file.path(OUTPUT_RESULTS_DIR, paste0("Differential_Expression_Results_", comp, ".csv"))
    write.csv(top_table_annotated, output_de_csv, row.names = FALSE)
    cat("  Differential expression results for", comp, "saved to:", output_de_csv, "\n")
  }
  cat("\nDifferential expression analysis completed for all comparisons.\n")
  
  # --- Step 6: Categorize Results for Post-DEA Enrichment Analysis ---
  cat("\n--- Categorizing proteins for Post-DEA enrichment analysis ---\n")
  significant_proteins_post_dea <- list()
  for (comp in comparisons) {
    df <- all_de_results[[comp]]
    significant_proteins_post_dea[[comp]]$over <- df %>% filter(adj.P.Val < 0.05 & logFC > 0) %>% pull(Protein_Accession)
    if (length(significant_proteins_post_dea[[comp]]$over) == 0) {
      cat(paste0("  No 'over-represented' proteins found for comparison ", comp, ".\n"))
    }
    significant_proteins_post_dea[[comp]]$under <- df %>% filter(adj.P.Val < 0.05 & logFC < 0) %>% pull(Protein_Accession)
    if (length(significant_proteins_post_dea[[comp]]$under) == 0) {
      cat(paste0("  No 'under-represented' proteins found for comparison ", comp, ".\n"))
    }
    cat(paste0("\nComparison ", comp, ":\n  Over-represented proteins: ", length(significant_proteins_post_dea[[comp]]$over), "\n  Under-represented proteins: ", length(significant_proteins_post_dea[[comp]]$under), "\n"))
  }
  
  # --- Step 7: Perform Post-DEA Enrichment Analysis (Modular Function Call) ---
  cat("\n--- Initiating Post-DEA enrichment analysis ---\n")
  
  for (comp in comparisons) {
    for (direction in c("OverRepresented", "UnderRepresented")) {
      gene_list <- if (direction == "OverRepresented") significant_proteins_post_dea[[comp]]$over else significant_proteins_post_dea[[comp]]$under
      
      # For GO Terms, iterate through categories
      go_categories <- c("biological_process", "molecular_function", "cellular_component")
      for (cat in go_categories) {
        perform_enrichment(
          gene_list = gene_list,
          background_genes = background_proteins,
          full_annotations_df = full_annotations_with_all_gos, # Use all GO rows for enrichment mapping
          interpro_raw_df = interpro_raw_df,
          annotation_type = "GO_Terms",
          comparison_name = comp,
          direction = direction,
          output_dir = OUTPUT_RESULTS_DIR,
          dpi = DPI,
          analysis_prefix = "Post_DEA",
          go_category = cat # Pass the specific GO category
        )
      }
      
      # For InterPro Methods, iterate through specific databases
      interpro_databases <- c("Pfam", "PANTHER", "SUPERFAMILY", "CDD", "Gene3D", "TIGRFAMs")
      for (db in interpro_databases) {
        perform_enrichment(
          gene_list = gene_list,
          background_genes = background_proteins,
          full_annotations_df = full_annotations_df_for_limma_and_background, # Still use distinct for this
          interpro_raw_df = interpro_raw_df, # Raw interpro data is crucial for granular enrichment
          annotation_type = "InterPro_Methods",
          comparison_name = comp,
          direction = direction,
          output_dir = OUTPUT_RESULTS_DIR,
          dpi = DPI,
          analysis_prefix = "Post_DEA",
          interpro_db = db # Pass the specific InterPro database
        )
      }
      
      # For EC Numbers and KO Numbers (general analysis)
      for (enrich_type in c("EC_Numbers", "KO_Numbers")) {
        perform_enrichment(
          gene_list = gene_list,
          background_genes = background_proteins,
          full_annotations_df = full_annotations_df_for_limma_and_background,
          interpro_raw_df = interpro_raw_df, # Keep as general parameter
          annotation_type = enrich_type,
          comparison_name = comp,
          direction = direction,
          output_dir = OUTPUT_RESULTS_DIR,
          dpi = DPI,
          analysis_prefix = "Post_DEA"
        )
      }
    }
  }
  cat("\n==============================================================================")
  cat("\nAll analyses (InterProScan, Data Consolidation, Pre-DEA FC-based Enrichment, limma, and Post-DEA Enrichment) completed.")
  cat("\nCheck the 'AnnotationPGresults_v3' folder for all output files.")
  cat("\n==============================================================================\n")
}

# Execute the main analysis workflow
main_analysis_workflow()