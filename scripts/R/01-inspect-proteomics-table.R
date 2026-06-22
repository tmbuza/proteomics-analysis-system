#!/usr/bin/env Rscript

# 01-inspect-proteomics-table.R
#
# Inspect a processed proteomics result table.
#
# Usage:
#   Rscript scripts/R/01-inspect-proteomics-table.R data/example/example-proteomics-results.csv results
#
# Optional:
#   Rscript scripts/R/01-inspect-proteomics-table.R data/input/proteomics-results.tsv results

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
      "Rscript scripts/R/01-inspect-proteomics-table.R <input_table> [output_dir]",
      sep = "\n"
    ),
    call. = FALSE
  )
}

input_file <- args[1]
output_dir <- ifelse(length(args) >= 2, args[2], "results")

if (!file.exists(input_file)) {
  stop(paste("Input file not found:", input_file), call. = FALSE)
}

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

summary_file <- file.path(output_dir, "table-inspection-summary.txt")
column_summary_file <- file.path(output_dir, "table-column-summary.tsv")
missing_values_file <- file.path(output_dir, "table-missing-values.tsv")
detected_columns_file <- file.path(output_dir, "table-detected-columns.tsv")

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

delimiter <- detect_delimiter(input_file)

proteomics_tbl <- readr::read_delim(
  file = input_file,
  delim = delimiter,
  show_col_types = FALSE,
  progress = FALSE,
  name_repair = "minimal"
)

column_names <- names(proteomics_tbl)

column_summary <- tibble(
  column_name = column_names,
  column_index = seq_along(column_names),
  column_type = vapply(proteomics_tbl, function(x) class(x)[1], character(1)),
  non_missing_values = vapply(proteomics_tbl, function(x) sum(!is.na(x) & x != ""), integer(1)),
  missing_values = vapply(proteomics_tbl, function(x) sum(is.na(x) | x == ""), integer(1)),
  missing_percent = round((missing_values / nrow(proteomics_tbl)) * 100, 2)
)

missing_values <- column_summary %>%
  select(column_name, missing_values, missing_percent) %>%
  arrange(desc(missing_values), column_name)

normalize_column <- function(x) {
  x %>%
    stringr::str_to_lower() %>%
    stringr::str_replace_all("[^a-z0-9]+", "_") %>%
    stringr::str_replace_all("^_|_$", "")
}

detect_columns <- function(column_names, patterns, category, exclude_patterns = character()) {
  normalized_columns <- normalize_column(column_names)
  combined_pattern <- paste(patterns, collapse = "|")

  detected <- tibble(
    category = category,
    column_name = column_names,
    normalized_column = normalized_columns,
    matched = stringr::str_detect(normalized_column, combined_pattern)
  )

  if (length(exclude_patterns) > 0) {
    exclude_pattern <- paste(exclude_patterns, collapse = "|")
    detected <- detected %>%
      mutate(excluded = stringr::str_detect(normalized_column, exclude_pattern)) %>%
      filter(matched, !excluded)
  } else {
    detected <- detected %>%
      filter(matched)
  }

  detected <- detected %>%
    select(category, column_name)

  if (nrow(detected) == 0) {
    detected <- tibble(
      category = category,
      column_name = "Not detected"
    )
  }

  detected
}

detected_columns <- bind_rows(
  detect_columns(
    column_names,
    c(
      "^protein_id$",
      "^protein_ids$",
      "^accession$",
      "^accessions$",
      "^uniprot$",
      "^uniprot_id$",
      "^uniprot_ids$",
      "^majority_protein_ids$",
      "^protein_group$",
      "^protein_groups$"
    ),
    "protein_identifier"
  ),
  detect_columns(
    column_names,
    c("^gene$", "^gene_name$", "^gene_symbol$", "^symbol$", "^genes$"),
    "gene_name"
  ),
  detect_columns(
    column_names,
    c("^protein_name$", "^protein_description$", "^description$", "^annotation$"),
    "protein_description"
  ),
  detect_columns(
    column_names,
    c("^log2fc$", "^log2_fc$", "^log2_fold_change$", "^logfc$", "^fold_change$", "^ratio$", "^difference$"),
    "fold_change"
  ),
  detect_columns(
    column_names,
    c("^p$", "^p_value$", "^pvalue$", "^p_val$", "^pval$", "^p_value_raw$", "^raw_p_value$"),
    "p_value",
    exclude_patterns = c("adjust", "^adj", "padj", "fdr", "q_value", "qvalue")
  ),
  detect_columns(
    column_names,
    c("^adjusted_p_value$", "^adj_p_value$", "^adj_p$", "^padj$", "^q_value$", "^qvalue$", "^fdr$", "^false_discovery_rate$"),
    "adjusted_p_value_or_fdr"
  ),
  detect_columns(
    column_names,
    c("intensity", "abundance", "lfq", "ibaq", "spectral_count", "area"),
    "abundance_or_intensity"
  ),
  detect_columns(
    column_names,
    c("peptide", "unique_peptide", "razor", "psm"),
    "peptide_or_evidence"
  )
)

readr::write_tsv(column_summary, column_summary_file)
readr::write_tsv(missing_values, missing_values_file)
readr::write_tsv(detected_columns, detected_columns_file)

detected_text <- detected_columns %>%
  group_by(category) %>%
  summarise(columns = paste(column_name, collapse = ", "), .groups = "drop") %>%
  mutate(line = paste0(category, ": ", columns)) %>%
  pull(line)

preview_text <- proteomics_tbl %>%
  head(6) %>%
  print(n = 6, width = Inf) %>%
  capture.output()

summary_lines <- c(
  "Proteomics Table Inspection Summary",
  "====================================",
  "",
  paste("Input file:", input_file),
  paste("Detected delimiter:", ifelse(delimiter == "\t", "tab", delimiter)),
  paste("Rows:", nrow(proteomics_tbl)),
  paste("Columns:", ncol(proteomics_tbl)),
  "",
  "Column names:",
  paste(column_names, collapse = ", "),
  "",
  "Detected candidate columns:",
  detected_text,
  "",
  "Output files:",
  paste("Column summary:", column_summary_file),
  paste("Missing values:", missing_values_file),
  paste("Detected columns:", detected_columns_file),
  "",
  "Preview:",
  preview_text
)

readr::write_lines(summary_lines, summary_file)

cat(paste(summary_lines, collapse = "\n"))
cat("\n")
cat("\nSummary written to:", summary_file, "\n")
