# Install and load necessary packages if they are not already installed
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
if (!requireNamespace("Biostrings", quietly = TRUE))
  BiocManager::install("Biostrings")
if (!requireNamespace("progress", quietly = TRUE))
  install.packages("progress")
if (!requireNamespace("dplyr", quietly = TRUE))
  install.packages("dplyr") # For data manipulation (e.g., rename, filter, select)
if (!requireNamespace("parallel", quietly = TRUE))
  install.packages("parallel") # For detecting CPU cores

library(Biostrings)
library(progress)
library(dplyr)
library(parallel) 

# --- Configuration Parameters ---

# 1. Path to your query FASTA file (amino acid sequences)
query_file <- "C:/Users/Alfonso/Desktop/AlfonsoOA_MSI/1Investigacion/5BiodegradacionPlasticos/11MetagenomaTenebrio/Functional_annotation/Functional_annotation/MAG_translated_genes/MetaBAT2Refined-0.2.faa"

# 2. Path to the BLASTP executable
blastp_path <- "C:/Users/Alfonso/Downloads/ncbi-blast-2.16.0+-x64-win64/ncbi-blast-2.16.0+/bin/blastp.exe" 

# 3. Path to the local BLAST database (prefix)
db_path <- "C:/Users/Alfonso/Desktop/AlfonsoOA_MSI/1Investigacion/5BiodegradacionPlasticos/10EnriquemientoKarol/Metagenomica/6bPlasticDBplus/plasticsdbplus"

# 4. Base directory for output
base_output_dir <- "C:/Users/Alfonso/Desktop/AlfonsoOA_MSI/1Investigacion/5BiodegradacionPlasticos/11MetagenomaTenebrio"
# Specific output folder name for BLAST results
output_folder_name <- "Tenebrio_Blastp"
output_folder <- file.path(base_output_dir, output_folder_name)

# 5. Dynamic number of threads for BLASTP
# Uses all but 1 core. You can change `max(1, num_available_cores - X)` where X is cores to reserve.
num_available_cores <- detectCores(logical = TRUE)
num_blast_threads <- max(1, num_available_cores - 1)
message(paste("Detected", num_available_cores, "CPU cores. Using", num_blast_threads, "threads for BLASTP."))

# 6. Default filtering thresholds for BLAST results
# IMPORTANT: Adjust these values here, and the filename will reflect them.
evalue_threshold_default <- 0.01
percent_identity_threshold_default <- 30
query_coverage_threshold_default <- 20

# 7. Experiment name for dynamic output filename
experiment_name <- "Tenebrio" # <--- IMPORTANT: Customize this name

# 8. Output filenames for raw results (always fixed)
raw_output_filename <- "raw_blastp_results.txt" # All raw results will be saved here with headers

# --- End of Configuration ---

# --- Pre-run Checks ---

message("\n--- Performing pre-run checks ---")

# Check if BLASTP executable exists
if (!file.exists(blastp_path)) {
  stop(paste("Error: BLASTP executable not found at", blastp_path, "\nPlease check the 'blastp_path' setting."))
} else {
  message(paste("BLASTP executable found:", blastp_path))
}

# Check if BLAST database files exist (checking for .phr, .pin, .psq as indicators)
db_files_exist <- all(file.exists(paste0(db_path, c(".phr", ".pin", ".psq"))))
if (!db_files_exist) {
  stop(paste("Error: BLAST database files not found for prefix", db_path, "\nEnsure the database has been properly formatted (e.g., with 'makeblastdb')."))
} else {
  message(paste("BLAST database files found for prefix:", db_path))
}

# Create the output folder if it doesn't exist
message(paste("Checking/creating output folder:", output_folder))
if (!dir.exists(output_folder)) {
  dir.create(output_folder, recursive = TRUE)
  message("Output folder created.")
} else {
  message("Output folder already exists.")
}

# Load query sequences
message(paste("Loading query sequences from:", query_file))
query_sequences <- readAAStringSet(query_file)
num_sequences <- length(query_sequences)
message(paste(num_sequences, "query sequences loaded."))

# --- BLAST Execution ---

# Path for the combined raw BLAST output
combined_raw_blast_output_path <- file.path(output_folder, raw_output_filename)

# Standard outfmt 6 column names for raw BLAST results
blast_outfmt6_cols <- c(
  "qseqid", "sseqid", "pident", "length", "mismatch", "gapopen",
  "qstart", "qend", "sstart", "send", "evalue", "bitscore", "stitle"
)

# Check if the combined raw BLAST output file already exists
# If it exists, skip the BLAST execution step to save time.
# This logic means if you want to force a re-run of BLAST, you must delete raw_blastp_results.txt
if (file.exists(combined_raw_blast_output_path) && file.info(combined_raw_blast_output_path)$size > 0) {
  message(paste("\nRaw BLAST results file found:", combined_raw_blast_output_path))
  message("Skipping BLAST execution and proceeding to result loading and filtering.")
  
  # Load the existing raw results with headers
  all_raw_blast_hits <- read.table(combined_raw_blast_output_path, header = TRUE, sep = "\t", stringsAsFactors = FALSE, quote = "") 
  
} else {
  message("\nRaw BLAST results file not found or empty. Starting BLASTP analysis...")
  message("This might take a while depending on the number of sequences and database size.")
  
  # List to temporarily store results for combination
  temp_raw_results_list <- list()
  
  # Progress bar for BLAST
  pb_blast <- progress_bar$new(
    format = "BLASTing [:bar] :percent eta: :eta",
    total = num_sequences,
    clear = FALSE
  )
  
  # Iterate over each query sequence for BLAST
  for (i in seq_along(query_sequences)) {
    # Create a temporary FASTA file for the current query sequence
    temp_fasta_file <- file.path(output_folder, paste0("temp_query_", i, ".fasta"))
    writeXStringSet(query_sequences[i], file = temp_fasta_file)
    
    # Command to execute BLASTP locally, outputting to stdout
    blast_command <- paste0(
      shQuote(blastp_path),
      " -query ", shQuote(temp_fasta_file),
      " -db ", shQuote(db_path),
      " -outfmt \"6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore stitle\"",
      " -num_threads ", num_blast_threads,
      " -comp_based_stats 2"
    )
    
    # Execute the BLAST command and capture its output
    blast_output <- system(blast_command, intern = TRUE, ignore.stderr = TRUE)
    
    if (length(blast_output) > 0) {
      # Convert output lines to a data frame
      current_results <- read.table(text = blast_output, header = FALSE, sep = "\t", stringsAsFactors = FALSE, quote = "")
      colnames(current_results) <- blast_outfmt6_cols # Assign column names
      temp_raw_results_list[[i]] <- current_results # Store using index for consistent order
    } else {
      # If no hits, store an empty data frame with correct column names to maintain structure
      temp_raw_results_list[[i]] <- data.frame(matrix(ncol = length(blast_outfmt6_cols), nrow = 0))
      colnames(temp_raw_results_list[[i]]) <- blast_outfmt6_cols
    }
    
    # Clean up temporary query file
    file.remove(temp_fasta_file)
    
    pb_blast$tick() # Advance progress bar
  }
  
  # Combine all raw results into a single data frame
  all_raw_blast_hits <- do.call(rbind, temp_raw_results_list)
  
  # Save the combined raw BLAST results to a file with headers
  if (!is.null(all_raw_blast_hits) && nrow(all_raw_blast_hits) > 0) {
    write.table(all_raw_blast_hits, file = combined_raw_blast_output_path, sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
    message(paste("\nAll raw BLAST results (with headers) saved to:", combined_raw_blast_output_path))
    message(paste("Total raw hits found:", nrow(all_raw_blast_hits)))
  } else {
    message("\nNo raw BLAST hits were found for any query sequence. The raw output file is empty.")
    # Create an empty file to indicate it was checked, and avoid re-running BLAST next time
    file.create(combined_raw_blast_output_path)
  }
}

message("\n--- Proceeding to filter BLAST results ---")

# --- Diagnóstico: Mostrar filtros actuales ---
message(paste("Filtros actuales: E-value <= ", evalue_threshold_default,
              ", Identidad >= ", percent_identity_threshold_default, "%",
              ", Cobertura >= ", query_coverage_threshold_default, "%"))

# Add query lengths to the raw hits for coverage calculation
# This step is crucial for calculating query coverage correctly
query_lengths_df <- data.frame(
  qseqid = names(query_sequences),
  query_length = width(query_sequences),
  stringsAsFactors = FALSE
)

# --- FIX CRÍTICO: Limpiar qseqid para asegurar que coincida con la salida de BLAST ---
# Elimina cualquier descripción después del primer espacio en los nombres de las queries FASTA
query_lengths_df$qseqid <- sub(" .*", "", query_lengths_df$qseqid)
message("Nombres de queries FASTA limpiados para coincidir con el formato qseqid de BLAST.")


# Ensure 'all_raw_blast_hits' is not empty before filtering
if (!is.null(all_raw_blast_hits) && nrow(all_raw_blast_hits) > 0) {
  filtered_blast_results <- all_raw_blast_hits %>%
    left_join(query_lengths_df, by = "qseqid") %>%
    # Calculate query coverage: ensure qstart/qend are ordered for length calculation
    mutate(
      qstart_adj = pmin(qstart, qend),
      qend_adj = pmax(qstart, qend),
      query_coverage = ((qend_adj - qstart_adj + 1) / query_length) * 100
    ) %>%
    # --- Diagnóstico: Chequeo de hits específicos ---
    {
      # Obtener las filas de los hits de interés antes de filtrar
      hits_to_check <- filter(., qseqid %in% c("OJEALBPC_00555", "OJEALBPC_00521"))
      if(nrow(hits_to_check) > 0) {
        message("\n--- Diagnóstico de hits específicos (ANTES de aplicar filtros): ---")
        for (i in 1:nrow(hits_to_check)) {
          hit <- hits_to_check[i, ]
          message(paste0("  Query: ", hit$qseqid, ", Subject: ", hit$sseqid))
          message(paste0("    E-value: ", hit$evalue, " (Umbral: ", evalue_threshold_default, ") -> Pasa: ", hit$evalue <= evalue_threshold_default))
          message(paste0("    Identidad: ", hit$pident, "% (Umbral: ", percent_identity_threshold_default, "%) -> Pasa: ", hit$pident >= percent_identity_threshold_default))
          message(paste0("    Longitud total query: ", hit$query_length, " (de FASTA)"))
          message(paste0("    Alineamiento qstart: ", hit$qstart, ", qend: ", hit$qend))
          message(paste0("    Cobertura calculada: ", hit$query_coverage, "% (Umbral: ", query_coverage_threshold_default, "%) -> Pasa: ", hit$query_coverage >= query_coverage_threshold_default))
          
          # Verificar si pasaría todos los filtros
          pass_all <- (hit$evalue <= evalue_threshold_default &&
                         hit$pident >= percent_identity_threshold_default &&
                         hit$query_coverage >= query_coverage_threshold_default)
          message(paste0("    ¿Pasa TODOS los filtros?: ", pass_all))
        }
        message("--- Fin del diagnóstico de hits específicos ---")
      } else {
        message("\nLos hits OJEALBPC_00555 y/o OJEALBPC_00521 no se encontraron en los resultados raw de BLAST.")
      }
      . # Asegura que el pipeline dplyr continúe con los datos originales
    } %>%
    # Apply filtering criteria using the default thresholds
    filter(
      evalue <= evalue_threshold_default,
      pident >= percent_identity_threshold_default,
      query_coverage >= query_coverage_threshold_default
    ) %>%
    # Select and reorder columns for the final output
    select(
      qseqid, sseqid, stitle, pident, length, mismatch, gapopen, qstart, qend, sstart, send,
      evalue, bitscore, query_length, query_coverage
    )
} else {
  message("Raw BLAST results table is empty. No filtering can be performed.")
  filtered_blast_results <- NULL
}

# --- IMPORTANT CHANGE HERE: Dynamically create filename just before saving ---
# Filtered filename includes current filter values and experiment name
# Use formatC for E-value to control scientific notation and decimal places
evalue_formatted <- formatC(evalue_threshold_default, format = "e", digits = 2) # e.g., 1.00e-02
# Adjust for common E-values like 0.01 to be displayed as 0.01, not 1e-02
if (evalue_threshold_default >= 0.001) { # For E-values like 0.01, 0.001 etc. display as decimals
  evalue_formatted <- formatC(evalue_threshold_default, format = "f", digits = max(2, -log10(evalue_threshold_default) + 1))
  evalue_formatted <- gsub("\\.", "p", evalue_formatted) # Replace . with p to be filename-safe
  evalue_formatted <- gsub("0+$", "", evalue_formatted) # Remove trailing zeros after decimal point
  evalue_formatted <- gsub("p$", "", evalue_formatted) # Remove trailing 'p' if no decimals
} else { # For very small E-values, keep scientific notation but make it filename-safe
  evalue_formatted <- gsub("e-0", "e-", evalue_formatted) # e.g., 1.00e-02 to 1.00e-2
  evalue_formatted <- gsub("\\.", "p", evalue_formatted)  # Replace . with p
  evalue_formatted <- gsub("-", "m", evalue_formatted)    # Replace - with m
}


filtered_output_filename_dynamic <- paste0(
  "E", evalue_formatted,
  "_I", percent_identity_threshold_default,
  "_C", query_coverage_threshold_default,
  "_blastp_", experiment_name, ".csv"
)
output_file_filtered <- file.path(output_folder, filtered_output_filename_dynamic)


if (!is.null(filtered_blast_results) && nrow(filtered_blast_results) > 0) {
  write.csv(filtered_blast_results, file = output_file_filtered, row.names = FALSE)
  message(paste("\nFiltered BLAST results saved to:", output_file_filtered))
  message(paste("  Using E-value <= ", evalue_threshold_default))
  message(paste("  Using Percent Identity >= ", percent_identity_threshold_default, "%"))
  message(paste("  Using Query Coverage >= ", query_coverage_threshold_default, "%"))
  message(paste("Total filtered hits found:", nrow(filtered_blast_results)))
} else {
  message("\nNo hits were found that met the specified filtering thresholds (",
          "E:", evalue_threshold_default, ", I:", percent_identity_threshold_default,
          ", C:", query_coverage_threshold_default, "%).")
  message("Consider adjusting 'evalue_threshold_default', 'percent_identity_threshold_default', or 'query_coverage_threshold_default' if you expect more distant matches.")
}

message("\n--- Script finished. ---")
message(paste("Raw BLAST results are available at:", combined_raw_blast_output_path))
message(paste("You can load this file manually in R (e.g., `read.table(\"", combined_raw_blast_output_path, "\", header = TRUE, sep = \"\\t\")`) and apply different filters without re-running BLAST."))