#!/usr/bin/env Rscript

# 10-visualize-proteomics-results.R
#
# Create core visualizations from differential proteomics results.
#
# Usage:
#   Rscript scripts/R/10-visualize-proteomics-results.R results/differential-proteins.tsv results 10

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(stringr)
  library(tibble)
  library(ggplot2)
})

args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 1) {
  stop(
    paste(
      "Usage:",
      "Rscript scripts/R/10-visualize-proteomics-results.R <differential_proteins.tsv> [output_dir] [top_n]",
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

figure_dir <- file.path(output_dir, "figures")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)

volcano_data_file <- file.path(output_dir, "volcano-plot-data.tsv")
top_plot_data_file <- file.path(output_dir, "top-significant-proteins-for-plot.tsv")
summary_file <- file.path(output_dir, "visualization-summary.tsv")

volcano_plot_file <- file.path(figure_dir, "volcano-plot.png")
top_plot_file <- file.path(figure_dir, "top-significant-proteins.png")
regulation_plot_file <- file.path(figure_dir, "regulation-summary.png")

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

diff_tbl <- readr::read_delim(
  file = input_file,
  delim = delimiter,
  show_col_types = FALSE,
  progress = FALSE,
  name_repair = "minimal"
)

column_names <- names(diff_tbl)

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

protein_label_col <- find_column(
  column_names,
  c("gene_symbol", "gene_symbol_clean", "gene_name", "protein_id", "primary_protein_id", "protein_name"),
  "protein label",
  required = FALSE
)

regulation_col <- find_column(
  column_names,
  c("regulation", "direction", "change_direction"),
  "regulation",
  required = FALSE
)

significance_col <- find_column(
  column_names,
  c("is_significant", "significant", "passes_filter"),
  "significance",
  required = FALSE
)

plot_tbl <- diff_tbl %>%
  mutate(
    log2fc_plot = as.numeric(.data[[log2fc_col]]),
    adjusted_p_value_plot = as.numeric(.data[[padj_col]]),
    adjusted_p_value_plot = ifelse(
      is.na(adjusted_p_value_plot) | adjusted_p_value_plot <= 0,
      NA_real_,
      adjusted_p_value_plot
    ),
    minus_log10_adjusted_p_value = -log10(adjusted_p_value_plot),
    protein_label = if (!is.na(protein_label_col)) {
      as.character(.data[[protein_label_col]])
    } else {
      paste0("protein_", row_number())
    },
    regulation_plot = if (!is.na(regulation_col)) {
      as.character(.data[[regulation_col]])
    } else {
      case_when(
        log2fc_plot > 0 ~ "upregulated",
        log2fc_plot < 0 ~ "downregulated",
        TRUE ~ "not_changed"
      )
    },
    is_significant_plot = if (!is.na(significance_col)) {
      as.logical(.data[[significance_col]])
    } else {
      !is.na(adjusted_p_value_plot)
    }
  ) %>%
  filter(!is.na(log2fc_plot), !is.na(minus_log10_adjusted_p_value))

top_tbl <- plot_tbl %>%
  filter(is_significant_plot) %>%
  mutate(abs_log2fc_plot = abs(log2fc_plot)) %>%
  arrange(adjusted_p_value_plot, desc(abs_log2fc_plot)) %>%
  slice_head(n = top_n)

regulation_summary <- plot_tbl %>%
  count(regulation_plot, name = "protein_count") %>%
  arrange(desc(protein_count))

volcano_plot <- ggplot(
  plot_tbl,
  aes(x = log2fc_plot, y = minus_log10_adjusted_p_value)
) +
  geom_point(aes(shape = regulation_plot)) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
  labs(
    title = "Volcano Plot",
    x = "log2 fold change",
    y = "-log10 adjusted p-value",
    shape = "Regulation"
  ) +
  theme_minimal(base_size = 12)

top_plot <- ggplot(
  top_tbl,
  aes(
    x = reorder(protein_label, log2fc_plot),
    y = log2fc_plot
  )
) +
  geom_col() +
  coord_flip() +
  labs(
    title = "Top Significant Proteins",
    x = "Protein",
    y = "log2 fold change"
  ) +
  theme_minimal(base_size = 12)

regulation_plot <- ggplot(
  regulation_summary,
  aes(x = regulation_plot, y = protein_count)
) +
  geom_col() +
  labs(
    title = "Regulation Summary",
    x = "Regulation class",
    y = "Protein count"
  ) +
  theme_minimal(base_size = 12)

ggsave(volcano_plot_file, volcano_plot, width = 7, height = 5, dpi = 300)
ggsave(top_plot_file, top_plot, width = 7, height = 5, dpi = 300)
ggsave(regulation_plot_file, regulation_plot, width = 6, height = 4, dpi = 300)

readr::write_tsv(plot_tbl, volcano_data_file)
readr::write_tsv(top_tbl, top_plot_data_file)

summary_tbl <- tibble(
  metric = c(
    "input_file",
    "input_rows",
    "plot_ready_rows",
    "top_n",
    "top_significant_proteins_plotted",
    "log2fc_column",
    "adjusted_p_value_column",
    "protein_label_column",
    "regulation_column"
  ),
  value = c(
    input_file,
    as.character(nrow(diff_tbl)),
    as.character(nrow(plot_tbl)),
    as.character(top_n),
    as.character(nrow(top_tbl)),
    log2fc_col,
    padj_col,
    ifelse(is.na(protein_label_col), "generated_labels", protein_label_col),
    ifelse(is.na(regulation_col), "derived_from_log2fc", regulation_col)
  )
)

readr::write_tsv(summary_tbl, summary_file)

cat("Proteomics visualization complete.\n")
cat("\n")
cat("Input file:", input_file, "\n")
cat("Plot-ready rows:", nrow(plot_tbl), "\n")
cat("Top proteins plotted:", nrow(top_tbl), "\n")
cat("\n")
cat("Output files:\n")
cat(" -", volcano_plot_file, "\n")
cat(" -", top_plot_file, "\n")
cat(" -", regulation_plot_file, "\n")
cat(" -", volcano_data_file, "\n")
cat(" -", top_plot_data_file, "\n")
cat(" -", summary_file, "\n")
