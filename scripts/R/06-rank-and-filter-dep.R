#!/usr/bin/env Rscript

# 06-rank-and-filter-dep.R
#
# Rank differentially abundant proteins from the filtered differential results table.
#
# Usage:
#   Rscript scripts/R/06-rank-and-filter-dep.R results/differential-proteins.tsv results 10

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
      "Rscript scripts/R/06-rank-and-filter-dep.R <differential_proteins.tsv> [output_dir] [top_n]",
      sep = "\n"
    ),
    call. = FALSE
  )
}

input_file <- args[1]
output_dir <- ifelse(length(args) >= 2, args[2], "results")
top_n <- ifelse(length(args) >= 3, as.integer(args[3]), 10)

if (!file.exists(input_file)) {
  stop(paste("Input file not found:", input_file), call. = FALSE)
}

if (is.na(top_n) || top_n < 1) {
  stop("top_n must be a positive integer.", call. = FALSE)
}

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

ranked_file <- file.path(output_dir, "ranked-proteins.tsv")
ranked_significant_file <- file.path(output_dir, "ranked-significant-proteins.tsv")
top_upregulated_file <- file.path(output_dir, "top-upregulated-proteins.tsv")
top_downregulated_file <- file.path(output_dir, "top-downregulated-proteins.tsv")
summary_file <- file.path(output_dir, "ranking-summary.tsv")

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

delimiter <- detect_delimiter(input_file)

dep_tbl <- readr::read_delim(
  file = input_file,
  delim = delimiter,
  show_col_types = FALSE,
  progress = FALSE,
  name_repair = "minimal"
)

column_names <- names(dep_tbl)

log2fc_col <- find_column(
  column_names,
  c("log2fc_numeric", "log2fc", "log2FC", "log2_fc", "log2_fold_change", "logFC"),
  "log2 fold change"
)

padj_col <- find_column(
  column_names,
  c("adjusted_p_value_numeric", "adjusted_p_value", "padj", "adj_p_value", "adj.P.Val", "q_value", "FDR"),
  "adjusted p-value"
)

abs_log2fc_col <- find_column(
  column_names,
  c("abs_log2fc", "absolute_log2fc", "abs_log2_fc"),
  "absolute log2 fold change",
  required = FALSE
)

is_significant_col <- find_column(
  column_names,
  c("is_significant", "significant", "passes_filter"),
  "significance flag",
  required = FALSE
)

regulation_col <- find_column(
  column_names,
  c("regulation", "direction", "change_direction"),
  "regulation direction",
  required = FALSE
)

ranked_tbl <- dep_tbl %>%
  mutate(
    .rank_log2fc = as.numeric(.data[[log2fc_col]]),
    .rank_padj = as.numeric(.data[[padj_col]]),
    .rank_abs_log2fc = if (!is.na(abs_log2fc_col)) {
      as.numeric(.data[[abs_log2fc_col]])
    } else {
      abs(.rank_log2fc)
    },
    .rank_is_significant = if (!is.na(is_significant_col)) {
      as.logical(.data[[is_significant_col]])
    } else {
      !is.na(.rank_padj)
    },
    .rank_regulation = if (!is.na(regulation_col)) {
      as.character(.data[[regulation_col]])
    } else {
      case_when(
        .rank_log2fc > 0 ~ "upregulated",
        .rank_log2fc < 0 ~ "downregulated",
        TRUE ~ "not_changed"
      )
    }
  ) %>%
  arrange(.rank_padj, desc(.rank_abs_log2fc)) %>%
  mutate(rank = row_number()) %>%
  relocate(rank)

ranked_significant_tbl <- ranked_tbl %>%
  filter(.rank_is_significant)

top_upregulated_tbl <- ranked_significant_tbl %>%
  filter(.rank_regulation == "upregulated") %>%
  arrange(.rank_padj, desc(.rank_abs_log2fc)) %>%
  slice_head(n = top_n)

top_downregulated_tbl <- ranked_significant_tbl %>%
  filter(.rank_regulation == "downregulated") %>%
  arrange(.rank_padj, desc(.rank_abs_log2fc)) %>%
  slice_head(n = top_n)

summary_tbl <- tibble(
  metric = c(
    "input_file",
    "total_proteins_ranked",
    "significant_proteins_ranked",
    "top_n_requested",
    "top_upregulated_exported",
    "top_downregulated_exported",
    "log2fc_column",
    "adjusted_p_value_column"
  ),
  value = c(
    input_file,
    as.character(nrow(ranked_tbl)),
    as.character(nrow(ranked_significant_tbl)),
    as.character(top_n),
    as.character(nrow(top_upregulated_tbl)),
    as.character(nrow(top_downregulated_tbl)),
    log2fc_col,
    padj_col
  )
)

ranked_tbl_export <- ranked_tbl %>%
  select(-starts_with(".rank_"))

ranked_significant_export <- ranked_significant_tbl %>%
  select(-starts_with(".rank_"))

top_upregulated_export <- top_upregulated_tbl %>%
  select(-starts_with(".rank_"))

top_downregulated_export <- top_downregulated_tbl %>%
  select(-starts_with(".rank_"))

readr::write_tsv(ranked_tbl_export, ranked_file)
readr::write_tsv(ranked_significant_export, ranked_significant_file)
readr::write_tsv(top_upregulated_export, top_upregulated_file)
readr::write_tsv(top_downregulated_export, top_downregulated_file)
readr::write_tsv(summary_tbl, summary_file)

cat("DEP ranking complete.\n")
cat("\n")
cat("Input file:", input_file, "\n")
cat("Total proteins ranked:", nrow(ranked_tbl), "\n")
cat("Significant proteins ranked:", nrow(ranked_significant_tbl), "\n")
cat("Top N:", top_n, "\n")
cat("\n")
cat("Output files:\n")
cat(" -", ranked_file, "\n")
cat(" -", ranked_significant_file, "\n")
cat(" -", top_upregulated_file, "\n")
cat(" -", top_downregulated_file, "\n")
cat(" -", summary_file, "\n")
