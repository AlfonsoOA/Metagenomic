# Metagenome Functional Data Analysis of *Tenebrio molitor*

This repository contains an R script (`metagenome_functional_analysis.R`) for the processing and statistical analysis of metagenome functional abundance data, specifically focusing on *Tenebrio molitor* samples. The script performs alpha and beta diversity analyses, variance partitioning, identification of differentially abundant genes for exploratory visualization, gene co-occurrence analysis, and hierarchical clustering.

## Project Description

The aim of this script is to provide a comprehensive workflow for analyzing gene abundance data (CPM) obtained from metagenomic analyses. It allows for the exploration of functional diversity patterns, differences between experimental conditions, and associations between genes.

## Directory Structure

The following directory structure is expected:
├── 0_Metagenome/
│   ├── 1_EC_relative_abundance_table_tenebrio_molitor.csv
│   └── Metagenome_Results/  (This directory will be created automatically for outputs)
├── metagenome_functional_analysis.R  (This script)
└── README.md
## Requirements

This script requires R (version 4.5.0 or higher) and the following R packages:

* `tidyverse`
* `vegan`
* `ggpubr`
* `reshape2`
* `pheatmap`
* `VennDiagram`
* `ape`
* `igraph`
* `psych`
* `cluster`
* `dendextend`
* `Hmisc`
* `microeco` (version 1.14.0 or higher)

You can install these packages by running the following command in your R console if you don't have them already:

```R
packages <- c("tidyverse", "vegan", "ggpubr", "reshape2", "pheatmap",
              "VennDiagram", "ape", "igraph", "psych", "cluster", "dendextend",
              "Hmisc", "microeco")
invisible(lapply(packages, function(x) if (!require(x, character.only = TRUE))) install.packages(x)))

Usage

    Place the script and data: Ensure that the 1_EC_relative_abundance_table_tenebrio_molitor.csv file is in the path specified by base_path within the script (e.g., C:/Users/Alfonso/Desktop/AlfonsoOA_MSI/1Investigacion/5BiodegradacionPlasticos/11MetagenomaTenebrio/0_Metagenome_Figure2).

    Configure base_path: Open the metagenome_functional_analysis.R script and adjust the base_path variable to the absolute path of your main project directory (0_Metagenome).
base_path <- "C:/path/to/your/0_Metagenome/directory" # Example path, replace with your actual path

Run the script: You can execute the entire script from RStudio or from the R command line:
source("metagenome_functional_analysis.R")
Generated Outputs

The script will create a Metagenome_Results directory within 0_Metagenome (or the path defined by base_path) and save the following files:

Tabular Data (.csv)

    alpha_diversity.csv: Table containing Shannon and Simpson diversity indices per sample.

    PCoA_coordinates.csv: Coordinates of samples in the PCoA space.

    PERMANOVA_results.csv: Detailed results of the PERMANOVA analysis.

    ANOSIM_results.csv: Results of the ANOSIM analysis.

    LEfSe_results.csv: Results from the LEfSe analysis, identifying differentially abundant biomarkers.

Figures (.png and .tiff)

    alpha_diversity_Shannon.png, alpha_diversity_Shannon.tiff: Boxplots of Shannon diversity by condition.

    alpha_diversity_Simpson.png, alpha_diversity_Simpson.tiff: Boxplots of Simpson diversity by condition.

    Rarefaction_curve.png, Rarefaction_curve.tiff: Rarefaction curves to assess sequencing depth.

    PCoA_plot.png, PCoA_plot.tiff: PCoA plot showing beta diversity, with PERMANOVA R² and p-value.

    NMDS_plot.png, NMDS_plot.tiff: NMDS plot showing beta diversity and model stress, with ANOSIM R and p-value.

    PERMANOVA_variance_pie_chart.png, PERMANOVA_variance_pie_chart.tiff: Pie chart visualizing the proportion of variance explained by the experimental condition and residual variance.

    heatmap_ANOVA_top_genes.png, heatmap_ANOVA_top_genes.tiff: Heatmap of the top 25 most abundant and ANOVA-significant genes.

    venn_diagram_genes.png, venn_diagram_genes.tiff: Venn diagram showing the overlap of genes present across different conditions.

    dendrogram_replicas.png, dendrogram_replicas.tiff: Dendrogram of hierarchical clustering of experimental replicates.

    gene_cooccurrence_network.png, gene_cooccurrence_network.tiff: Gene co-occurrence network based on Spearman correlation.

    LEfSe_biomarkers.png, LEfSe_biomarkers.tiff: Bar plot visualizing differentially abundant functional biomarkers identified by LEfSe.

Analyses Performed

The script executes the following analyses in the specified order:

    Data Loading and Preparation: Imports abundance data and defines the condition factor.

    Alpha Diversity: Calculates Shannon and Simpson indices, and generates boxplots.

    Sample Rarefaction: Generates rarefaction curves based on rounded CPM data, utilizing the minimum total counts as sampling depth.

    Beta Diversity (Ordination and Statistical Tests):

        PCoA (Principal Coordinate Analysis): Explores structural differences in gene abundance profiles. Calculates and reports the percentage of variance explained by the first two principal coordinates.

        PERMANOVA (Permutational Multivariate Analysis of Variance): Statistically tests for significant differences in gene abundance profiles among experimental conditions using Bray-Curtis dissimilarity and 999 permutations. Quantifies the proportion of total variance explained (R2) by the "Condition" factor.

        NMDS (Non-metric Multidimensional Scaling): Explores structural differences in gene abundance profiles. Performs a two-dimensional ordination and reports the stress value.

        ANOSIM (Analysis of Similarities): Statistically tests for significant differences in gene abundance profiles among experimental conditions.

    Variance Partitioning (PERMANOVA): Visualizes the proportion of variance explained by the "Condition" factor versus residual variance using a pie chart.

    Hierarchical Clustering: Applies hierarchical clustering to samples based on Bray-Curtis dissimilarity using Ward's method (ward.D2).

    Exploratory Heatmap of Differentially Abundant Genes (ANOVA): Performs one-way ANOVA per gene, adjusts p-values (FDR), selects the top 25 most abundant significant genes (adjusted $p \< 0.05$), and visualizes their abundance patterns in a row-scaled heatmap with hierarchical clustering of samples and genes.

    Venn Diagram: Constructs a Venn diagram to assess the overlap of prevalent genes, considering a gene "present" if its sum abundance across all replicates for a condition exceeds 0 CPM.

    Gene Co-occurrence Network Analysis: Constructs a network of statistically significant positive and negative associations between gene abundances using Spearman's rank correlation (absolute correlation geq0.7, adjusted p-value leq0.01). Visualizes the network with igraph and a Fruchterman-Reingold layout.

    LEfSe (Linear Discriminant Analysis Effect Size) Analysis: Identifies differentially abundant functional biomarkers using trans_diff from the microeco package, with FDR correction applied to LEfSe p-values.

Contact

For any questions or feedback, please contact Alfonso Olaya Abril (b22olaba@uco.es).
