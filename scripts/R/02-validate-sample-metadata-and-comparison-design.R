#!/usr/bin/env Rscript

# 03-validate-sample-metadata-and-comparison-design.R
#
# Validate sample metadata and comparison design for a results-first proteomics workflow.
#
# Usage:
#   Rscript scripts/R/03-validate-sample-metadata-and-comparison-design.R \
#     data/example/example-proteomics-results.csv \
#     data/example/example-sample-metadata.csv \
#     data/example/example-comparison-design.csv \
#     results

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(stringr)
  library(tibble)
})

args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 3) {
  stop(
    paste(
      "Usage:",
      "Rscript scripts/R/03-validate-sample-metadata-and-comparison-design.R <proteomics_results> <sample_metadata> <comparison_design> [output_dir]",
      sep = "\n"
    ),
    call. = FALSE
  )
}

results_file <- args[1]
metadata_file <- args[2]
comparison_file <- args[3]
output_dir <- ifelse(length(args) >= 4, args[4], "results")

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

validation_summary_file <- file.path(output_dir, "metadata-validation-summary.tsv")
condition_counts_file <- file.path(output_dir, "condition-sample-counts.tsv")
comparison_validated_file <- file.path(output_dir, "comparison-design-validated.tsv")

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

read_table_auto <- function(file) {
  readr::read_delim(
    file = file,
    delim = detect_delimiter(file),
    show_col_types = FALSE,
    progress = FALSE,
    name_repair = "minimal"
  )
}

normalize_column <- function(x) {
  x %>%
    stringr::str_to_lower() %>%
    stringr::str_replace_all("[^a-z0-9]+", "_") %>%
    stringr::str_replace_all("^_|_$", "")
}

has_required_columns <- function(tbl, required_columns) {
  all(normalize_column(required_columns) %in% normalize_column(names(tbl)))
}

get_column_name <- function(tbl, target_column) {
  normalized_names <- normalize_column(names(tbl))
  normalized_target <- normalize_column(target_column)

  idx <- match(normalized_target, normalized_names)

  if (is.na(idx)) {
    return(NA_character_)
  }

  names(tbl)[idx]
}

validation_summary <- tibble(
  check = character(),
  status = character(),
  detail = character()
)

add_check <- function(summary_tbl, check, passed, detail) {
  bind_rows(
    summary_tbl,
    tibble(
      check = check,
      status = ifelse(passed, "PASS", "FAIL"),
      detail = detail
    )
  )
}

results_exists <- file.exists(results_file)
metadata_exists <- file.exists(metadata_file)
comparison_exists <- file.exists(comparison_file)

validation_summary <- add_check(validation_summary, "proteomics_results_file_exists", results_exists, results_file)
validation_summary <- add_check(validation_summary, "sample_metadata_file_exists", metadata_exists, metadata_file)
validation_summary <- add_check(validation_summary, "comparison_design_file_exists", comparison_exists, comparison_file)

if (!results_exists || !metadata_exists || !comparison_exists) {
  readr::write_tsv(validation_summary, validation_summary_file)
  stop("One or more required input files are missing. See metadata-validation-summary.tsv.", call. = FALSE)
}

results_tbl <- read_table_auto(results_file)
metadata_tbl <- read_table_auto(metadata_file)
comparison_tbl <- read_table_auto(comparison_file)

required_metadata_columns <- c("sample_id", "condition")
required_comparison_columns <- c(
  "comparison",
  "numerator_condition",
  "denominator_condition",
  "positive_log2fc_interpretation"
)

metadata_required_present <- has_required_columns(metadata_tbl, required_metadata_columns)
comparison_required_present <- has_required_columns(comparison_tbl, required_comparison_columns)

validation_summary <- add_check(
  validation_summary,
  "required_sample_metadata_columns_present",
  metadata_required_present,
  paste(required_metadata_columns, collapse = ", ")
)

validation_summary <- add_check(
  validation_summary,
  "required_comparison_design_columns_present",
  comparison_required_present,
  paste(required_comparison_columns, collapse = ", ")
)

if (!metadata_required_present || !comparison_required_present) {
  readr::write_tsv(validation_summary, validation_summary_file)
  stop("Required metadata or comparison-design columns are missing. See metadata-validation-summary.tsv.", call. = FALSE)
}

condition_col <- get_column_name(metadata_tbl, "condition")
comparison_col <- get_column_name(comparison_tbl, "comparison")
numerator_col <- get_column_name(comparison_tbl, "numerator_condition")
denominator_col <- get_column_name(comparison_tbl, "denominator_condition")

condition_counts <- metadata_tbl %>%
  count(.data[[condition_col]], name = "sample_count") %>%
  rename(condition = 1) %>%
  arrange(condition)

metadata_conditions <- unique(metadata_tbl[[condition_col]])

comparison_validated <- comparison_tbl %>%
  mutate(
    numerator_found_in_metadata = .data[[numerator_col]] %in% metadata_conditions,
    denominator_found_in_metadata = .data[[denominator_col]] %in% metadata_conditions,
    comparison_conditions_valid = numerator_found_in_metadata & denominator_found_in_metadata
  )

all_comparison_conditions_valid <- all(comparison_validated$comparison_conditions_valid)

validation_summary <- add_check(
  validation_summary,
  "comparison_conditions_found_in_metadata",
  all_comparison_conditions_valid,
  paste("Metadata conditions:", paste(metadata_conditions, collapse = ", "))
)

if ("comparison" %in% normalize_column(names(results_tbl))) {
  results_comparison_col <- names(results_tbl)[match("comparison", normalize_column(names(results_tbl)))]
  result_comparisons <- unique(results_tbl[[results_comparison_col]])
  design_comparisons <- unique(comparison_tbl[[comparison_col]])

  comparisons_match <- all(result_comparisons %in% design_comparisons)

  validation_summary <- add_check(
    validation_summary,
    "comparison_labels_match_result_table",
    comparisons_match,
    paste(
      "Result comparisons:",
      paste(result_comparisons, collapse = ", "),
      "| Design comparisons:",
      paste(design_comparisons, collapse = ", ")
    )
  )
} else {
  validation_summary <- add_check(
    validation_summary,
    "comparison_labels_match_result_table",
    FALSE,
    "No comparison column detected in proteomics results table."
  )
}

readr::write_tsv(validation_summary, validation_summary_file)
readr::write_tsv(condition_counts, condition_counts_file)
readr::write_tsv(comparison_validated, comparison_validated_file)

cat("Sample metadata and comparison design validation complete.\n")
cat("\n")
cat("Proteomics results:", results_file, "\n")
cat("Sample metadata:", metadata_file, "\n")
cat("Comparison design:", comparison_file, "\n")
cat("\n")
cat("Output files:\n")
cat(" -", validation_summary_file, "\n")
cat(" -", condition_counts_file, "\n")
cat(" -", comparison_validated_file, "\n")
