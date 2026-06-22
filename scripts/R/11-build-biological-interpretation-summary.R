#!/usr/bin/env Rscript

# 11-build-biological-interpretation-summary.R
#
# Build a structured biological interpretation summary from proteomics workflow outputs.
#
# Usage:
#   Rscript scripts/R/11-build-biological-interpretation-summary.R results results/biological-interpretation-summary.md

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(stringr)
  library(tibble)
})

args <- commandArgs(trailingOnly = TRUE)

results_dir <- ifelse(length(args) >= 1, args[1], "results")
output_file <- ifelse(
  length(args) >= 2,
  args[2],
  file.path(results_dir, "biological-interpretation-summary.md")
)

if (!dir.exists(results_dir)) {
  stop(paste("Results directory not found:", results_dir), call. = FALSE)
}

dir.create(dirname(output_file), recursive = TRUE, showWarnings = FALSE)

read_tsv_if_exists <- function(path) {
  if (file.exists(path)) {
    return(readr::read_tsv(path, show_col_types = FALSE, progress = FALSE))
  }

  tibble()
}

metric_lookup <- function(tbl, metric_name, default = "not_available") {
  if (!all(c("metric", "value") %in% names(tbl))) {
    return(default)
  }

  value <- tbl %>%
    filter(.data$metric == metric_name) %>%
    pull(.data$value)

  if (length(value) == 0) {
    return(default)
  }

  as.character(value[1])
}

differential_summary <- read_tsv_if_exists(file.path(results_dir, "differential-summary.tsv"))
ranking_summary <- read_tsv_if_exists(file.path(results_dir, "ranking-summary.tsv"))
identifier_summary <- read_tsv_if_exists(file.path(results_dir, "identifier-cleaning-summary.tsv"))
enrichment_summary <- read_tsv_if_exists(file.path(results_dir, "enrichment-input-summary.tsv"))
string_summary <- read_tsv_if_exists(file.path(results_dir, "string-network-summary.tsv"))
top_proteins <- read_tsv_if_exists(file.path(results_dir, "top-significant-proteins-for-plot.tsv"))
ranked_significant <- read_tsv_if_exists(file.path(results_dir, "ranked-significant-proteins.tsv"))

figure_dir <- file.path(results_dir, "figures")
figure_files <- c(
  "volcano-plot.png",
  "top-significant-proteins.png",
  "regulation-summary.png"
)

available_figures <- figure_files[file.exists(file.path(figure_dir, figure_files))]

format_table_preview <- function(tbl, max_rows = 10) {
  if (nrow(tbl) == 0) {
    return(c("No table available."))
  }

  preview <- tbl %>% slice_head(n = max_rows)
  capture.output(print(preview, n = max_rows, width = Inf))
}

differential_lines <- c()

if (nrow(differential_summary) > 0) {
  differential_lines <- c(
    differential_lines,
    "The differential abundance summary file was detected.",
    "",
    "```text",
    capture.output(print(differential_summary, n = Inf, width = Inf)),
    "```"
  )
} else {
  differential_lines <- c(
    differential_lines,
    "Differential summary file was not detected."
  )
}

top_protein_lines <- c()

if (nrow(top_proteins) > 0) {
  likely_label_columns <- intersect(
    c("protein_label", "gene_symbol", "gene_symbol_clean", "protein_id", "primary_protein_id", "protein_name"),
    names(top_proteins)
  )

  label_col <- ifelse(length(likely_label_columns) > 0, likely_label_columns[1], names(top_proteins)[1])

  top_names <- top_proteins %>%
    slice_head(n = 10) %>%
    pull(.data[[label_col]]) %>%
    as.character()

  top_protein_lines <- c(
    top_protein_lines,
    paste("Top proteins available for review:", paste(top_names, collapse = ", "))
  )
} else if (nrow(ranked_significant) > 0) {
  top_protein_lines <- c(
    top_protein_lines,
    "Ranked significant protein table was detected. Review the top rows for biological interpretation."
  )
} else {
  top_protein_lines <- c(
    top_protein_lines,
    "No ranked significant protein table was detected."
  )
}

identifier_lines <- c(
  paste("Input rows:", metric_lookup(identifier_summary, "input_rows")),
  paste("Unique cleaned protein identifiers:", metric_lookup(identifier_summary, "unique_cleaned_protein_identifiers")),
  paste("Gene symbols available:", metric_lookup(identifier_summary, "gene_symbols_available")),
  paste("Multi-identifier rows:", metric_lookup(identifier_summary, "multi_identifier_rows")),
  paste("Possible contaminant rows:", metric_lookup(identifier_summary, "possible_contaminant_rows")),
  paste("Possible reverse or decoy rows:", metric_lookup(identifier_summary, "possible_reverse_or_decoy_rows"))
)

enrichment_lines <- c(
  paste("All gene identifiers:", metric_lookup(enrichment_summary, "all_gene_identifiers")),
  paste("Upregulated gene identifiers:", metric_lookup(enrichment_summary, "upregulated_gene_identifiers")),
  paste("Downregulated gene identifiers:", metric_lookup(enrichment_summary, "downregulated_gene_identifiers")),
  paste("All protein identifiers:", metric_lookup(enrichment_summary, "all_protein_identifiers")),
  paste("Upregulated protein identifiers:", metric_lookup(enrichment_summary, "upregulated_protein_identifiers")),
  paste("Downregulated protein identifiers:", metric_lookup(enrichment_summary, "downregulated_protein_identifiers"))
)

string_lines <- c(
  paste("Usable unique STRING identifiers:", metric_lookup(string_summary, "usable_unique_identifiers")),
  paste("Upregulated identifiers:", metric_lookup(string_summary, "upregulated_identifiers")),
  paste("Downregulated identifiers:", metric_lookup(string_summary, "downregulated_identifiers")),
  paste("STRING identifier column:", metric_lookup(string_summary, "string_identifier_column"))
)

figure_lines <- if (length(available_figures) > 0) {
  paste("- `", file.path("figures", available_figures), "`", sep = "")
} else {
  "No figure files were detected."
}

summary_lines <- c(
  "# Biological Interpretation Summary",
  "",
  "## Analysis Overview",
  "",
  "This summary was generated from the available proteomics workflow outputs.",
  "",
  paste("Results directory:", results_dir),
  "",
  "## Differential Abundance Summary",
  "",
  differential_lines,
  "",
  "## Top-Ranked Proteins",
  "",
  top_protein_lines,
  "",
  "## Identifier Readiness",
  "",
  paste("- ", identifier_lines, sep = ""),
  "",
  "## Functional Enrichment Readiness",
  "",
  paste("- ", enrichment_lines, sep = ""),
  "",
  "## STRING Network Readiness",
  "",
  paste("- ", string_lines, sep = ""),
  "",
  "## Available Figures",
  "",
  figure_lines,
  "",
  "## Interpretation Notes",
  "",
  "Use the ranked proteins, enrichment-ready identifiers, STRING-ready identifiers, and figures to write a biological interpretation.",
  "",
  "Recommended interpretation pattern:",
  "",
  "```text",
  "The differential abundance results show ...",
  "The top-ranked proteins include ...",
  "The enrichment-ready gene/protein lists suggest follow-up analysis of ...",
  "The STRING-ready identifiers can be used to test whether these proteins form a connected network ...",
  "```",
  "",
  "## Limitations",
  "",
  "- Interpretation depends on the quality of the upstream proteomics result table.",
  "- Missing or ambiguous identifiers may reduce annotation and enrichment coverage.",
  "- Protein groups and isoforms may require manual review.",
  "- Functional enrichment and STRING results should be interpreted with organism and database context.",
  "- Differential abundance does not prove mechanism without additional evidence.",
  "",
  "## Next Steps",
  "",
  "- Review top-ranked proteins manually.",
  "- Run organism-specific GO/pathway enrichment if not already completed.",
  "- Run STRING network analysis using the prepared upload lists.",
  "- Compare enrichment and network results for coherent biological themes.",
  "- Integrate figures and interpretation into the final reproducible report."
)

readr::write_lines(summary_lines, output_file)

cat("Biological interpretation summary created.\n")
cat("\n")
cat("Output file:", output_file, "\n")
