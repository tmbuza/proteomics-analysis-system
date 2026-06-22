#!/usr/bin/env Rscript

# 04-qc-result-table.R
#
# Perform quality control checks on a processed proteomics result table.
#
# Usage:
#   Rscript scripts/R/04-qc-result-table.R data/example/example-proteomics-results.csv results

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
      "Rscript scripts/R/04-qc-result-table.R <input_table> [output_dir]",
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

qc_summary_file <- file.path(output_dir, "result-table-qc-summary.tsv")
column_qc_file <- file.path(output_dir, "result-table-column-qc.tsv")
duplicate_file <- file.path(output_dir, "result-table-duplicate-proteins.tsv")
problem_rows_file <- file.path(output_dir, "result-table-problem-rows.tsv")
qc_report_file <- file.path(output_dir, "result-table-qc-report.txt")

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

get_column_name <- function(tbl, candidate_names) {
  normalized_names <- normalize_column(names(tbl))
  normalized_candidates <- normalize_column(candidate_names)

  matched_index <- match(normalized_candidates, normalized_names)
  matched_index <- matched_index[!is.na(matched_index)]

  if (length(matched_index) == 0) {
    return(NA_character_)
  }

  names(tbl)[matched_index[1]]
}

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

delimiter <- detect_delimiter(input_file)

result_tbl <- readr::read_delim(
  file = input_file,
  delim = delimiter,
  show_col_types = FALSE,
  progress = FALSE,
  name_repair = "minimal"
)

protein_id_col <- get_column_name(result_tbl, c("protein_id", "Protein_ID", "accession", "Accession", "uniprot_id", "UniProt_ID"))
log2fc_col <- get_column_name(result_tbl, c("log2fc", "log2FC", "log2_fc", "log2_fold_change", "logFC"))
p_value_col <- get_column_name(result_tbl, c("p_value", "pvalue", "p_val", "pval", "P.Value"))
padj_col <- get_column_name(result_tbl, c("adjusted_p_value", "padj", "adj_p_value", "adj.P.Val", "q_value", "FDR"))
comparison_col <- get_column_name(result_tbl, c("comparison", "contrast", "condition_comparison"))

required_columns <- tibble(
  logical_name = c("protein_id", "log2fc", "p_value", "adjusted_p_value", "comparison"),
  detected_column = c(protein_id_col, log2fc_col, p_value_col, padj_col, comparison_col)
)

required_columns_present <- all(!is.na(required_columns$detected_column))

column_qc <- tibble(
  column_name = names(result_tbl),
  column_type = vapply(result_tbl, function(x) class(x)[1], character(1)),
  missing_values = vapply(result_tbl, function(x) sum(is.na(x) | x == ""), integer(1)),
  missing_percent = round((missing_values / nrow(result_tbl)) * 100, 2),
  unique_values = vapply(result_tbl, function(x) length(unique(x)), integer(1))
)

qc_summary <- tibble(
  check = character(),
  status = character(),
  detail = character()
)

qc_summary <- add_check(qc_summary, "table_readable", TRUE, input_file)
qc_summary <- add_check(qc_summary, "rows_present", nrow(result_tbl) > 0, paste("Rows:", nrow(result_tbl)))
qc_summary <- add_check(qc_summary, "columns_present", ncol(result_tbl) > 0, paste("Columns:", ncol(result_tbl)))
qc_summary <- add_check(
  qc_summary,
  "required_columns_present",
  required_columns_present,
  paste(required_columns$logical_name, required_columns$detected_column, sep = "=", collapse = "; ")
)

protein_missing <- NA_integer_
duplicates_tbl <- tibble()

if (!is.na(protein_id_col)) {
  protein_values <- result_tbl[[protein_id_col]]
  protein_missing <- sum(is.na(protein_values) | protein_values == "")

  duplicates_tbl <- result_tbl %>%
    filter(!is.na(.data[[protein_id_col]]) & .data[[protein_id_col]] != "") %>%
    count(.data[[protein_id_col]], name = "duplicate_count") %>%
    filter(duplicate_count > 1) %>%
    rename(protein_id = 1) %>%
    arrange(desc(duplicate_count), protein_id)

  qc_summary <- add_check(
    qc_summary,
    "protein_identifiers_present",
    protein_missing == 0,
    paste("Missing protein IDs:", protein_missing)
  )

  qc_summary <- add_check(
    qc_summary,
    "duplicate_protein_identifiers_absent",
    nrow(duplicates_tbl) == 0,
    paste("Duplicated protein IDs:", nrow(duplicates_tbl))
  )
} else {
  duplicates_tbl <- tibble(protein_id = character(), duplicate_count = integer())
  qc_summary <- add_check(qc_summary, "protein_identifiers_present", FALSE, "Protein ID column not detected.")
  qc_summary <- add_check(qc_summary, "duplicate_protein_identifiers_absent", FALSE, "Protein ID column not detected.")
}

check_numeric_column <- function(tbl, column_name) {
  if (is.na(column_name)) {
    return(FALSE)
  }

  values <- suppressWarnings(as.numeric(tbl[[column_name]]))
  all(!is.na(values) | is.na(tbl[[column_name]]) | tbl[[column_name]] == "")
}

log2fc_numeric <- check_numeric_column(result_tbl, log2fc_col)
p_value_numeric <- check_numeric_column(result_tbl, p_value_col)
padj_numeric <- check_numeric_column(result_tbl, padj_col)

qc_summary <- add_check(qc_summary, "log2fc_numeric", log2fc_numeric, paste("Column:", log2fc_col))
qc_summary <- add_check(qc_summary, "p_value_numeric", p_value_numeric, paste("Column:", p_value_col))
qc_summary <- add_check(qc_summary, "adjusted_p_value_numeric", padj_numeric, paste("Column:", padj_col))

comparison_missing <- NA_integer_

if (!is.na(comparison_col)) {
  comparison_missing <- sum(is.na(result_tbl[[comparison_col]]) | result_tbl[[comparison_col]] == "")
  qc_summary <- add_check(
    qc_summary,
    "comparison_labels_present",
    comparison_missing == 0,
    paste("Missing comparison labels:", comparison_missing)
  )
} else {
  qc_summary <- add_check(qc_summary, "comparison_labels_present", FALSE, "Comparison column not detected.")
}

problem_rows <- result_tbl %>%
  mutate(.row_number = row_number())

if (!is.na(protein_id_col)) {
  problem_rows <- problem_rows %>%
    mutate(.missing_protein_id = is.na(.data[[protein_id_col]]) | .data[[protein_id_col]] == "")
} else {
  problem_rows <- problem_rows %>%
    mutate(.missing_protein_id = TRUE)
}

if (!is.na(log2fc_col)) {
  problem_rows <- problem_rows %>%
    mutate(.missing_log2fc = is.na(.data[[log2fc_col]]) | .data[[log2fc_col]] == "")
} else {
  problem_rows <- problem_rows %>%
    mutate(.missing_log2fc = TRUE)
}

if (!is.na(padj_col)) {
  problem_rows <- problem_rows %>%
    mutate(.missing_adjusted_p_value = is.na(.data[[padj_col]]) | .data[[padj_col]] == "")
} else {
  problem_rows <- problem_rows %>%
    mutate(.missing_adjusted_p_value = TRUE)
}

if (!is.na(comparison_col)) {
  problem_rows <- problem_rows %>%
    mutate(.missing_comparison = is.na(.data[[comparison_col]]) | .data[[comparison_col]] == "")
} else {
  problem_rows <- problem_rows %>%
    mutate(.missing_comparison = TRUE)
}

problem_rows <- problem_rows %>%
  filter(
    .missing_protein_id |
      .missing_log2fc |
      .missing_adjusted_p_value |
      .missing_comparison
  )

ready_for_filtering <- all(qc_summary$status == "PASS")

qc_summary <- add_check(
  qc_summary,
  "ready_for_filtering",
  ready_for_filtering,
  ifelse(ready_for_filtering, "All required QC checks passed.", "One or more QC checks failed.")
)

report_lines <- c(
  "Proteomics Result Table QC Report",
  "=================================",
  "",
  paste("Input file:", input_file),
  paste("Detected delimiter:", ifelse(delimiter == "\t", "tab", delimiter)),
  paste("Rows:", nrow(result_tbl)),
  paste("Columns:", ncol(result_tbl)),
  "",
  "Detected key columns:",
  paste("protein_id:", protein_id_col),
  paste("log2fc:", log2fc_col),
  paste("p_value:", p_value_col),
  paste("adjusted_p_value:", padj_col),
  paste("comparison:", comparison_col),
  "",
  "QC checks:",
  paste(qc_summary$check, qc_summary$status, qc_summary$detail, sep = " | "),
  "",
  "Output files:",
  paste("QC summary:", qc_summary_file),
  paste("Column QC:", column_qc_file),
  paste("Duplicate proteins:", duplicate_file),
  paste("Problem rows:", problem_rows_file)
)

readr::write_tsv(qc_summary, qc_summary_file)
readr::write_tsv(column_qc, column_qc_file)
readr::write_tsv(duplicates_tbl, duplicate_file)
readr::write_tsv(problem_rows, problem_rows_file)
readr::write_lines(report_lines, qc_report_file)

cat(paste(report_lines, collapse = "\n"))
cat("\n")
