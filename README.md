# Metagenomic
Metagenomic and ASV Analysis of Tenebrio molitor Gut Microbiome
This repository contains scripts and associated documentation for the bioinformatic analysis of 16S rRNA gene (ASVs) and whole metagenome sequencing data. The aim of this study is to investigate methodological limitations in identifying plastic degradation genes in the Tenebrio molitor gut microbiome.
Repository Structure
├── README.md ├── scripts/ │ ├── script1.R # 16S Microbial Diversity Statistical Analysis │ ├── script2.sh # InterProScan Pipeline for Functional Annotations │ ├── script3.R # Metagenome Gene Abundance Data Processing and Statistical Analysis │ ├── script4.R # Abundance Visualization and Enrichment Analysis │ └── script5.R # BLAST Results Filtering ├── data/ │ ├── (raw_sequences/) # Directory for raw sequences (not included, accessible via NCBI SRA) │ └── (processed_data/) # Directory for intermediate processed data ├── results/ │ └── figures/ │ ├── Figure1_ASVs_2872025.jpg │ ├── Figure2_Metagenome_28072025.jpg │ ├── Figure3_AnnotationsAccumulatedCPM_28072025.jpg │ ├── FigureS2_AccumulatedPlotsTop15.jpg │ ├── FigureS3_RarefractionCurves.jpg │ └── FigureS4_PostDEA-Dotplots.jpg │ └── tables/ # Directory for result tables
Requirements
This project requires the following tools and libraries. The use of a virtual environment (such as Conda) for dependency management is highly recommended.
•	R (version 4.5.0 or higher):
o	vegan
o	ggplot2
o	patchwork
o	pairwiseAdonis
o	indicspecies
o	pheatmap
o	VennDiagram
o	igraph
o	microeco (version 1.14.0)
o	limma
o	clusterProfiler
•	Bash:
o	InterProScan (version 5.75-106.0 or higher)
o	BLAST+ (version 2.16.0+ or higher)
•	Other Bioinformatics Software:
o	Qiime2 (version 2024.2) (for 16S rRNA gene amplicon sequencing bioinformatics)
o	Cutadapt
o	FastQC
o	DADA2
o	PICRUSt2
o	STAMP (version 2.1.3)
o	Bowtie2
o	Kraken 2
o	KronaTools
o	nf-core/mag pipeline (version 3.2.1)
o	Fastp
o	MEGAHIT
o	QUAST
o	MetaBAT2
o	MaxBin2
o	DAS Tool
o	CheckM
o	CAT/BAT (v5.2.3)
o	Prodigal
o	DIAMOND
o	Salmon
o	Prokka
o	OmicBox
o	GhostKOALA
Script Usage
Below is a description of each script included in this repository, forming part of the metagenome and 16S ASV data analysis workflow.
1. script1.R: 16S Microbial Community Diversity Statistical Analysis
•	Purpose: This R script performs multivariate statistical analyses of microbial community diversity from 16S rRNA gene sequencing data. It calculates community ecology metrics, assesses alpha and beta diversity, and conducts indicator species analysis. It also infers predicted functional profiles using PICRUSt2 and analyzes differentially abundant pathways with STAMP.
•	Technology: R (version 4.5.0)
•	Key Dependencies: vegan, indicspecies, pairwiseAdonis.
•	Input:
o	ASV (Amplicon Sequence Variant) count table.
o	ASV taxonomic assignment data (used by PICRUSt2).
•	Output:
o	Alpha and beta diversity results.
o	PERMANOVA and ISA results.
o	Predicted pathway abundance tables.
o	Related statistical plots.
2. script2.sh: InterProScan Functional Annotation Pipeline
•	Purpose: This Bash script implements an automated pipeline for performing protein domain and functional site annotations using InterProScan. It handles verification, installation (if necessary), parallel execution of InterProScan, and post-processing of results to generate comprehensive annotation tables.
•	Technology: Bash Script
•	Key Dependencies: InterProScan.
•	Usage: The script likely takes a FASTA protein file as input and may require parameters for the number of CPU cores to use and the InterProScan applications to enable.
•	Input: FASTA files of protein sequences.
•	Output:
o	Comprehensive annotation tables, including domain presence/absence matrices.
o	GO terms and pathway associations for proteins.
3. script3.R: Metagenome Gene Abundance Data Processing and Statistical Analysis
•	Purpose: This R script processes and analyzes gene abundance data (CPM) derived from metagenomic sequencing. It calculates alpha diversity metrics, generates rarefaction curves, performs ordinations (NMDS, PCoA), conducts statistical tests (PERMANOVA, ANOSIM), hierarchical clustering, heatmap visualization using one-way ANOVA, Venn diagram construction, gene co-occurrence network analysis, and LEfSe analysis.
•	Technology: R (version 4.5.0)
•	Key Dependencies: microeco (for core analyses), ggplot2, VennDiagram, igraph.
•	Input: Gene abundance table (CPM normalized).
•	Output:
o	Alpha diversity results and rarefaction curves.
o	Ordination plots (NMDS, PCoA).
o	PERMANOVA and ANOSIM test results.
o	Heatmaps and Venn diagrams.
o	Gene co-occurrence networks.
o	LEfSe results.
4. script4.R: Abundance Visualization and Enrichment Analysis
•	Purpose: This R script focuses on visualizing abundance data and performing functional enrichment analyses from diverse annotation datasets consolidated into a single dataframe. It includes pre-differential expression analysis (Pre-DEA), differential expression analysis (DEA) using the limma package, and functional enrichment analysis (both Pre-DEA and Post-DEA) using clusterProfiler.
•	Technology: R
•	Key Dependencies: limma, clusterProfiler.
•	Input: Consolidated dataframe of protein annotations with abundance data.
•	Output:
o	Abundance plots and visualizations.
o	Differential expression analysis results.
o	Functional enrichment results (e.g., enriched GO terms, pathways).
5. script5.R: BLAST Results Filtering
•	Purpose: This R script is designed to filter raw BLAST output results against a custom database. It applies stringent criteria such as E-value, percentage identity, and query coverage to refine the results and obtain highly relevant matches.
•	Technology: R
•	Key Dependencies: Standard R data manipulation packages.
•	Input: Raw BLAST output file (e.g., in tabular format).
•	Output: Filtered BLAST results table.
Key Results and Visualizations
Some of the key results generated by the scripts in this project are visually represented by the following figures:
•	ASV Community Analysis:
•	Metagenome Community Structure:
•	Accumulated CPM Annotations:
•	Top 15 Accumulated Plots:
•	Rarefaction Curves:
•	Post-DEA Dotplots:
References
This work makes use of the following tools and methodologies:
•	R Core Team. (2024). R: A language and environment for statistical computing. R Foundation for Statistical Computing, Vienna, Austria.
•	Callahan, B. J., McMurdie, P. J., Rosen, M. J., Han, A. W., Johnson, A. J. A., & Holmes, S. P. (2016). DADA2: High-resolution sample inference from Illumina amplicon data. Nature Methods, 13(7), 581-583.
•	InterProScan: A sequence analysis application that combines different protein signature recognition methods into one resource.
•	Wilcoxon, F., & Wilcox, R. A. (1964). Some rapid approximate statistical procedures. Lederle Laboratories.
•	Parks, D. H., Imelfort, M., Skennerton, C. T., Hugenholtz, P., & Tyson, G. W. (2015). CheckM: assessing the quality of microbial genomes recovered from isolates, single cells, and metagenomes. Genome Research, 25(7), 1043-1055.
•	Seemann, T. (2014). Prokka: rapid prokaryotic genome annotation. Bioinformatics, 30(14), 2068-2069.
•	Buchfink, B., Xie, C., & Huson, D. H. (2015). Fast and sensitive protein alignment using DIAMOND. Nature Methods, 12(1), 59-60.
•	Langmead, B., & Salzberg, S. L. (2012). Fast gapped-read alignment with Bowtie2. Nature Methods, 9(4), 357-359.
•	Bray, J. R., & Curtis, J. T. (1957). An ordination of the upland forest communities of southern Wisconsin. Ecological Monographs, 27(4), 325-349.
•	Anderson, M. J. (2001). A new method for non-parametric multivariate analysis of variance. Austral Ecology, 26(1), 32-46.
•	De Cáceres, M., & Legendre, P. (2009). Indicspecies: a R package for the analysis of species associations. Journal of Vegetation Science, 20(3), 760-765.
•	Douglas, G. M., Maffei, V. J., Zaneveld, J. R., Yurgel, S. N., Brown, J. R., Taylor, C. M., ... & Huttenhower, C. (2019). PICRUSt2 for prediction of metagenome functions. Nature Biotechnology, 37(12), 1459-1464.
•	Parks, D. H., Tyson, G. W., Hugenholtz, P., & Beiko, R. G. (2014). STAMP: statistical analysis of taxonomic and functional profiles. Bioinformatics, 30(20), 3122-3124.
•	Wang, Q., Garrity, G. M., Tiedje, J. M., & Cole, J. R. (2007). Naive Bayesian classifier for rapid assignment of rRNA sequences into new or existing taxa. Applied and Environmental Microbiology, 73(16), 5261-5267.
•	White, J. R., Nagarajan, N., & Pop, M. (2009). Statistical approaches for detecting differentially abundant features in clinical metagenomic samples. PLoS Computational Biology, 5(7), e1000372.
For a complete list of references, please refer to the References.docx file.
Contact
For any questions or comments, please contact Alfonso Olaya Abril (b22olaba@uco.es,  AlfonsoOA).
