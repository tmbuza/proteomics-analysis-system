#!/usr/bin/env Rscript

# 09-prepare-string-network-input.R
#
# Prepare STRING network input tables and upload lists from cleaned proteomics identifiers.
#
# Usage:
#   Rscript scripts/R/09-prepare-string-network-input.R results/string-network-input.tsv results

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
      "Rscript scripts/R/09-prepare-string-network-input.R <string_network_input.tsv> [output_dir]",
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

string_all_file <- file.path(output_dir, "string-input-all.tsv")
string_up_file <- file.path(output_dir, "string-input-upregulated.tsv")
string_down_file <- file.path(output_dir, "string-input-downregulated.tsv")

upload_all_file <- file.path(output_dir, "string-upload-list-all.txt")
upload_up_file <- file.path(output_dir, "string-upload-list-upregulated.txt")
upload_down_file <- file.path(output_dir, "string-upload-list-downregulated.txt")

summary_file <- file.path(output_dir, "string-network-summary.tsv")

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

find_column <- function(column_names, candidate_names, label, required = TRUE) {
  normalized_columns <- normalize_column(column_names)
  normalized_candidates <- normalize_column(candidate_names)

  matched_index <- match(normalized_candidates, normalized_columns)
  matched_index <- matched_index[!is.na(matched_index)]

  if (length(matched_index) == 0) {
    if (required) {
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
    } else {
      return(NA_character_)
    }
  }

  column_names[matched_index[1]]
}

write_upload_list <- function(tbl, id_col, output_file) {
  ids <- tbl %>%
    transmute(identifier = stringr::str_trim(as.character(.data[[id_col]]))) %>%
    filter(!is.na(identifier), identifier != "") %>%
    distinct(identifier) %>%
    arrange(identifier) %>%
    pull(identifier)

  readr::write_lines(ids, output_file)
  length(ids)
}

delimiter <- detect_delimiter(input_file)

string_tbl <- readr::read_delim(
  file = input_file,
  delim = delimiter,
  show_col_types = FALSE,
  progress = FALSE,
  name_repair = "minimal"
)

column_names <- names(string_tbl)

string_id_col <- find_column(
  column_names,
  c("string_identifier", "gene_symbol_clean", "gene_symbol", "primary_protein_id", "protein_id"),
  "STRING identifier"
)

regulation_col <- find_column(
  column_names,
  c("regulation", "direction", "change_direction"),
  "regulation direction",
  required = FALSE
)

rank_col <- find_column(
  column_names,
  c("rank", "protein_rank"),
  "rank",
  required = FALSE
)

prepared_tbl <- string_tbl %>%
  mutate(
    string_identifier_clean = stringr::str_trim(as.character(.data[[string_id_col]]))
  ) %>%
  filter(!is.na(string_identifier_clean), string_identifier_clean != "") %>%
  distinct(string_identifier_clean, .keep_all = TRUE)

if (!is.na(rank_col)) {
  prepared_tbl <- prepared_tbl %>%
    arrange(.data[[rank_col]])
} else {
  prepared_tbl <- prepared_tbl %>%
    arrange(string_identifier_clean)
}

if (!is.na(regulation_col)) {
  up_tbl <- prepared_tbl %>%
    filter(.data[[regulation_col]] == "upregulated")

  down_tbl <- prepared_tbl %>%
    filter(.data[[regulation_col]] == "downregulated")
} else {
  up_tbl <- prepared_tbl[0, , drop = FALSE]
  down_tbl <- prepared_tbl[0, , drop = FALSE]
}

all_upload_count <- write_upload_list(prepared_tbl, "string_identifier_clean", upload_all_file)
up_upload_count <- write_upload_list(up_tbl, "string_identifier_clean", upload_up_file)
down_upload_count <- write_upload_list(down_tbl, "string_identifier_clean", upload_down_file)

readr::write_tsv(prepared_tbl, string_all_file)
readr::write_tsv(up_tbl, string_up_file)
readr::write_tsv(down_tbl, string_down_file)

summary_tbl <- tibble(
  metric = c(
    "input_file",
    "input_rows",
    "string_identifier_column",
    "regulation_column",
    "usable_unique_identifiers",
    "upregulated_identifiers",
    "downregulated_identifiers",
    "missing_or_empty_identifiers_removed",
    "upload_list_all_count",
    "upload_list_upregulated_count",
    "upload_list_downregulated_count"
  ),
  value = c(
    input_file,
    as.character(nrow(string_tbl)),
    string_id_col,
    ifelse(is.na(regulation_col), "not_detected", regulation_col),
    as.character(nrow(prepared_tbl)),
    as.character(nrow(up_tbl)),
    as.character(nrow(down_tbl)),
    as.character(nrow(string_tbl) - nrow(prepared_tbl)),
    as.character(all_upload_count),
    as.character(up_upload_count),
    as.character(down_upload_count)
  )
)

readr::write_tsv(summary_tbl, summary_file)

cat("STRING network input preparation complete.\n")
cat("\n")
cat("Input file:", input_file, "\n")
cat("Input rows:", nrow(string_tbl), "\n")
cat("Usable unique identifiers:", nrow(prepared_tbl), "\n")
cat("STRING identifier column:", string_id_col, "\n")
cat("\n")
cat("Output files:\n")
cat(" -", string_all_file, "\n")
cat(" -", string_up_file, "\n")
cat(" -", string_down_file, "\n")
cat(" -", upload_all_file, "\n")
cat(" -", upload_up_file, "\n")
cat(" -", upload_down_file, "\n")
cat(" -", summary_file, "\n")
