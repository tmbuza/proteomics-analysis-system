#!/usr/bin/env Rscript

# 02-filter-differential-proteins.R
#
# Filter differential proteins from a processed proteomics result table.
#
# Usage:
#   Rscript scripts/R/02-filter-differential-proteins.R data/example/example-proteomics-results.csv results 1 0.05
#
# Arguments:
#   1. input table
#   2. output directory
#   3. absolute log2 fold-change threshold
#   4. adjusted p-value threshold

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(stringr)
  library(tibble)
})

args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 1) {
  stop(
    paste(
      "Usage:",
      "Rscript scripts/R/02-filter-differential-proteins.R <input_table> [output_dir] [abs_log2fc_threshold] [adjusted_p_value_threshold]",
      sep = "\n"
    ),
    call. = FALSE
  )
}

input_file <- args[1]
output_dir <- ifelse(length(args) >= 2, args[2], "results")
abs_log2fc_threshold <- ifelse(length(args) >= 3, as.numeric(args[3]), 1)
adjusted_p_value_threshold <- ifelse(length(args) >= 4, as.numeric(args[4]), 0.05)

if (!file.exists(input_file)) {
  stop(paste("Input file not found:", input_file), call. = FALSE)
}

if (is.na(abs_log2fc_threshold)) {
  stop("The absolute log2 fold-change threshold must be numeric.", call. = FALSE)
}

if (is.na(adjusted_p_value_threshold)) {
  stop("The adjusted p-value threshold must be numeric.", call. = FALSE)
}

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

differential_file <- file.path(output_dir, "differential-proteins.tsv")
significant_file <- file.path(output_dir, "significant-proteins.tsv")
upregulated_file <- file.path(output_dir, "upregulated-proteins.tsv")
downregulated_file <- file.path(output_dir, "downregulated-proteins.tsv")
summary_file <- file.path(output_dir, "differential-summary.tsv")

detect_delimiter <- function(file) {
  first_line <- readr::read_lines(file, n_max = 1)

  if (stringr::str_detect(first_line, "\t")) {
    return("\t")
  }

  if (stringr::str_detect(first_line, ",")) {
    return(",")
  }

  return("\t")
}

normalize_column <- function(x) {
  x %>%
    stringr::str_to_lower() %>%
    stringr::str_replace_all("[^a-z0-9]+", "_") %>%
    stringr::str_replace_all("^_|_$", "")
}

find_required_column <- function(column_names, candidate_names, label) {
  normalized_columns <- normalize_column(column_names)
  normalized_candidates <- normalize_column(candidate_names)

  matched_index <- match(normalized_candidates, normalized_columns)
  matched_index <- matched_index[!is.na(matched_index)]

  if (length(matched_index) == 0) {
    stop(
      paste(
        "Required column not found for:", label,
        "\nAccepted names:",
        paste(candidate_names, collapse = ", "),
        "\nAvailable columns:",
        paste(column_names, collapse = ", ")
      ),
      call. = FALSE
    )
  }

  column_names[matched_index[1]]
}

delimiter <- detect_delimiter(input_file)

proteomics_tbl <- readr::read_delim(
  file = input_file,
  delim = delimiter,
  show_col_types = FALSE,
  progress = FALSE,
  name_repair = "minimal"
)

column_names <- names(proteomics_tbl)

log2fc_col <- find_required_column(
  column_names,
  c("log2fc", "log2FC", "log2_fc", "log2_fold_change", "logFC"),
  "log2 fold change"
)

padj_col <- find_required_column(
  column_names,
  c("adjusted_p_value", "padj", "adj_p_value", "adj.P.Val", "q_value", "FDR"),
  "adjusted p-value or FDR"
)

protein_id_col <- find_required_column(
  column_names,
  c("protein_id", "Protein_ID", "accession", "Accession", "uniprot_id", "UniProt_ID"),
  "protein identifier"
)

comparison_col <- NULL
comparison_candidates <- c("comparison", "contrast", "condition_comparison")
comparison_matches <- match(normalize_column(comparison_candidates), normalize_column(column_names))
comparison_matches <- comparison_matches[!is.na(comparison_matches)]

if (length(comparison_matches) > 0) {
  comparison_col <- column_names[comparison_matches[1]]
}

differential_tbl <- proteomics_tbl %>%
  mutate(
    log2fc_numeric = as.numeric(.data[[log2fc_col]]),
    adjusted_p_value_numeric = as.numeric(.data[[padj_col]]),
    abs_log2fc = abs(log2fc_numeric),
    passes_padj = adjusted_p_value_numeric <= adjusted_p_value_threshold,
    passes_log2fc = abs_log2fc >= abs_log2fc_threshold,
    is_significant = passes_padj & passes_log2fc,
    regulation = case_when(
      is_significant & log2fc_numeric > 0 ~ "upregulated",
      is_significant & log2fc_numeric < 0 ~ "downregulated",
      TRUE ~ "not_significant"
    )
  ) %>%
  arrange(adjusted_p_value_numeric, desc(abs_log2fc))

significant_tbl <- differential_tbl %>%
  filter(is_significant)

upregulated_tbl <- differential_tbl %>%
  filter(regulation == "upregulated")

downregulated_tbl <- differential_tbl %>%
  filter(regulation == "downregulated")

if (!is.null(comparison_col)) {
  summary_tbl <- differential_tbl %>%
    group_by(comparison = .data[[comparison_col]]) %>%
    summarise(
      total_proteins = n(),
      significant_proteins = sum(is_significant, na.rm = TRUE),
      upregulated_proteins = sum(regulation == "upregulated", na.rm = TRUE),
      downregulated_proteins = sum(regulation == "downregulated", na.rm = TRUE),
      not_significant_proteins = sum(regulation == "not_significant", na.rm = TRUE),
      abs_log2fc_threshold = abs_log2fc_threshold,
      adjusted_p_value_threshold = adjusted_p_value_threshold,
      .groups = "drop"
    )
} else {
  summary_tbl <- tibble(
    comparison = "not_provided",
    total_proteins = nrow(differential_tbl),
    significant_proteins = sum(differential_tbl$is_significant, na.rm = TRUE),
    upregulated_proteins = sum(differential_tbl$regulation == "upregulated", na.rm = TRUE),
    downregulated_proteins = sum(differential_tbl$regulation == "downregulated", na.rm = TRUE),
    not_significant_proteins = sum(differential_tbl$regulation == "not_significant", na.rm = TRUE),
    abs_log2fc_threshold = abs_log2fc_threshold,
    adjusted_p_value_threshold = adjusted_p_value_threshold
  )
}

readr::write_tsv(differential_tbl, differential_file)
readr::write_tsv(significant_tbl, significant_file)
readr::write_tsv(upregulated_tbl, upregulated_file)
readr::write_tsv(downregulated_tbl, downregulated_file)
readr::write_tsv(summary_tbl, summary_file)

cat("Differential abundance filtering complete.\n")
cat("\n")
cat("Input file:", input_file, "\n")
cat("Rows:", nrow(proteomics_tbl), "\n")
cat("Protein ID column:", protein_id_col, "\n")
cat("Log2FC column:", log2fc_col, "\n")
cat("Adjusted p-value column:", padj_col, "\n")
cat("Absolute log2FC threshold:", abs_log2fc_threshold, "\n")
cat("Adjusted p-value threshold:", adjusted_p_value_threshold, "\n")
cat("\n")
cat("Output files:\n")
cat(" -", differential_file, "\n")
cat(" -", significant_file, "\n")
cat(" -", upregulated_file, "\n")
cat(" -", downregulated_file, "\n")
cat(" -", summary_file, "\n")
