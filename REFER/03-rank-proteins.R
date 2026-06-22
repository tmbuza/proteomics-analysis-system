#!/usr/bin/env Rscript

# 03-rank-proteins.R
#
# Rank proteins from a differential proteomics result table.
#
# Usage:
#   Rscript scripts/R/03-rank-proteins.R results/differential-proteins.tsv results 10

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(tibble)
})

args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 1) {
  stop(
    paste(
      "Usage:",
      "Rscript scripts/R/03-rank-proteins.R <differential_table> [output_dir] [top_n]",
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

proteins_tbl <- readr::read_tsv(
  input_file,
  show_col_types = FALSE,
  progress = FALSE
)

required_columns <- c(
  "log2fc_numeric",
  "adjusted_p_value_numeric",
  "abs_log2fc",
  "is_significant",
  "regulation"
)

missing_columns <- setdiff(required_columns, names(proteins_tbl))

if (length(missing_columns) > 0) {
  stop(
    paste(
      "Required columns are missing from the differential table:",
      paste(missing_columns, collapse = ", "),
      "\nRun scripts/R/02-filter-differential-proteins.R first."
    ),
    call. = FALSE
  )
}

ranked_tbl <- proteins_tbl %>%
  mutate(
    significance_rank_group = if_else(is_significant, 1L, 2L)
  ) %>%
  arrange(
    significance_rank_group,
    adjusted_p_value_numeric,
    desc(abs_log2fc)
  ) %>%
  mutate(rank_overall = row_number()) %>%
  select(rank_overall, everything(), -significance_rank_group)

ranked_significant_tbl <- ranked_tbl %>%
  filter(is_significant) %>%
  arrange(adjusted_p_value_numeric, desc(abs_log2fc)) %>%
  mutate(rank_significant = row_number()) %>%
  select(rank_significant, everything())

top_upregulated_tbl <- ranked_tbl %>%
  filter(regulation == "upregulated") %>%
  arrange(adjusted_p_value_numeric, desc(log2fc_numeric)) %>%
  slice_head(n = top_n) %>%
  mutate(rank_upregulated = row_number()) %>%
  select(rank_upregulated, everything())

top_downregulated_tbl <- ranked_tbl %>%
  filter(regulation == "downregulated") %>%
  arrange(adjusted_p_value_numeric, log2fc_numeric) %>%
  slice_head(n = top_n) %>%
  mutate(rank_downregulated = row_number()) %>%
  select(rank_downregulated, everything())

summary_tbl <- tibble(
  total_proteins = nrow(ranked_tbl),
  significant_proteins = nrow(ranked_significant_tbl),
  upregulated_proteins = sum(ranked_tbl$regulation == "upregulated", na.rm = TRUE),
  downregulated_proteins = sum(ranked_tbl$regulation == "downregulated", na.rm = TRUE),
  not_significant_proteins = sum(ranked_tbl$regulation == "not_significant", na.rm = TRUE),
  top_n_requested = top_n,
  top_upregulated_exported = nrow(top_upregulated_tbl),
  top_downregulated_exported = nrow(top_downregulated_tbl)
)

readr::write_tsv(ranked_tbl, ranked_file)
readr::write_tsv(ranked_significant_tbl, ranked_significant_file)
readr::write_tsv(top_upregulated_tbl, top_upregulated_file)
readr::write_tsv(top_downregulated_tbl, top_downregulated_file)
readr::write_tsv(summary_tbl, summary_file)

cat("Protein ranking complete.\n\n")
cat("Input file:", input_file, "\n")
cat("Rows ranked:", nrow(ranked_tbl), "\n")
cat("Significant proteins:", nrow(ranked_significant_tbl), "\n")
cat("Top N requested:", top_n, "\n\n")
cat("Output files:\n")
cat(" -", ranked_file, "\n")
cat(" -", ranked_significant_file, "\n")
cat(" -", top_upregulated_file, "\n")
cat(" -", top_downregulated_file, "\n")
cat(" -", summary_file, "\n")
