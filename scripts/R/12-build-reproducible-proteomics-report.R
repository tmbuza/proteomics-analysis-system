#!/usr/bin/env Rscript

# 12-build-reproducible-proteomics-report.R
#
# Build and render a reproducible proteomics report from workflow outputs.
#
# Usage:
#   Rscript scripts/R/12-build-reproducible-proteomics-report.R results reports

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(stringr)
  library(tibble)
})

args <- commandArgs(trailingOnly = TRUE)

results_dir <- ifelse(length(args) >= 1, args[1], "results")
report_dir <- ifelse(length(args) >= 2, args[2], "reports")

if (!dir.exists(results_dir)) {
  stop(paste("Results directory not found:", results_dir), call. = FALSE)
}

dir.create(report_dir, recursive = TRUE, showWarnings = FALSE)

report_qmd <- file.path(report_dir, "proteomics-report.qmd")
report_html <- file.path(report_dir, "proteomics-report.html")

read_tsv_if_exists <- function(path) {
  if (file.exists(path)) {
    return(readr::read_tsv(path, show_col_types = FALSE, progress = FALSE))
  }

  tibble()
}

file_status_line <- function(path) {
  if (file.exists(path)) {
    paste("- Found:", path)
  } else {
    paste("- Missing:", path)
  }
}

table_block <- function(title, path, max_rows = 10) {
  if (!file.exists(path)) {
    return(c(
      paste0("## ", title),
      "",
      paste("File not available:", path),
      ""
    ))
  }

  tbl <- read_tsv_if_exists(path)

  if (nrow(tbl) == 0) {
    return(c(
      paste0("## ", title),
      "",
      paste("File available but contains no rows:", path),
      ""
    ))
  }

  preview_path <- path
  printed <- capture.output(print(tbl %>% slice_head(n = max_rows), n = max_rows, width = Inf))

  c(
    paste0("## ", title),
    "",
    paste("Source file:", preview_path),
    "",
    "```text",
    printed,
    "```",
    ""
  )
}

image_block <- function(title, path, alt_text) {
  if (!file.exists(path)) {
    return(c(
      paste0("## ", title),
      "",
      paste("Figure not available:", path),
      ""
    ))
  }

  rel_path <- file.path("..", path)

  c(
    paste0("## ", title),
    "",
    paste0("![", alt_text, "](", rel_path, ")"),
    ""
  )
}

interpretation_path <- file.path(results_dir, "biological-interpretation-summary.md")

interpretation_lines <- if (file.exists(interpretation_path)) {
  readr::read_lines(interpretation_path)
} else {
  c(
    "## Biological Interpretation",
    "",
    paste("Biological interpretation summary was not found:", interpretation_path),
    "",
    "Run:",
    "",
    "```bash",
    "Rscript scripts/R/11-build-biological-interpretation-summary.R results results/biological-interpretation-summary.md",
    "```"
  )
}

input_files <- c(
  file.path(results_dir, "differential-summary.tsv"),
  file.path(results_dir, "ranked-significant-proteins.tsv"),
  file.path(results_dir, "identifier-cleaning-summary.tsv"),
  file.path(results_dir, "enrichment-input-summary.tsv"),
  file.path(results_dir, "string-network-summary.tsv"),
  file.path(results_dir, "visualization-summary.tsv"),
  interpretation_path,
  file.path(results_dir, "figures", "volcano-plot.png"),
  file.path(results_dir, "figures", "top-significant-proteins.png"),
  file.path(results_dir, "figures", "regulation-summary.png")
)

session_lines <- capture.output(sessionInfo())

report_lines <- c(
  "---",
  "title: "Proteomics Analysis Report"",
  "subtitle: "Results-First Differential Protein Abundance Workflow"",
  "format:",
  "  html:",
  "    toc: true",
  "    toc-depth: 3",
  "    number-sections: true",
  "    code-fold: true",
  "execute:",
  "  echo: false",
  "  warning: false",
  "  message: false",
  "---",
  "",
  "# Executive Summary",
  "",
  "This report summarizes a results-first proteomics analysis workflow.",
  "",
  "The report was generated from structured workflow outputs including differential abundance summaries, ranked protein tables, identifier-cleaning outputs, enrichment-ready inputs, STRING-ready inputs, figures, and a biological interpretation summary.",
  "",
  "# Input and Output Overview",
  "",
  paste("Report generated:", as.character(Sys.time())),
  "",
  paste("Results directory:", results_dir),
  "",
  paste("Report directory:", report_dir),
  "",
  "## Input File Status",
  "",
  vapply(input_files, file_status_line, character(1)),
  "",
  table_block("Differential Abundance Summary", file.path(results_dir, "differential-summary.tsv")),
  table_block("Ranked Significant Proteins", file.path(results_dir, "ranked-significant-proteins.tsv")),
  table_block("Identifier Cleaning Summary", file.path(results_dir, "identifier-cleaning-summary.tsv")),
  table_block("Functional Enrichment Readiness", file.path(results_dir, "enrichment-input-summary.tsv")),
  table_block("STRING Network Readiness", file.path(results_dir, "string-network-summary.tsv")),
  table_block("Visualization Summary", file.path(results_dir, "visualization-summary.tsv")),
  "# Figures",
  "",
  image_block(
    "Volcano Plot",
    file.path(results_dir, "figures", "volcano-plot.png"),
    "Volcano plot"
  ),
  image_block(
    "Top Significant Proteins",
    file.path(results_dir, "figures", "top-significant-proteins.png"),
    "Top significant proteins"
  ),
  image_block(
    "Regulation Summary",
    file.path(results_dir, "figures", "regulation-summary.png"),
    "Regulation summary"
  ),
  "# Biological Interpretation",
  "",
  interpretation_lines,
  "",
  "# Limitations",
  "",
  "- This report is based on processed proteomics result tables.",
  "- Raw mass spectrometry files were not reprocessed in this results-first edition.",
  "- Interpretation depends on the quality of the upstream differential abundance analysis.",
  "- Identifier mapping and annotation coverage can affect enrichment and network analysis.",
  "- STRING and GO interpretation require organism-specific context.",
  "- Differential protein abundance does not prove mechanism without additional validation.",
  "",
  "# Reproducibility Notes",
  "",
  "The report was generated from workflow outputs saved in the results directory.",
  "",
  "## R Session Information",
  "",
  "```text",
  session_lines,
  "```",
  ""
)

readr::write_lines(report_lines, report_qmd)

render_status <- "not_attempted"

if (requireNamespace("quarto", quietly = TRUE)) {
  render_status <- tryCatch(
    {
      quarto::quarto_render(
        input = report_qmd,
        output_file = basename(report_html),
        quiet = TRUE
      )
      "rendered_with_quarto_r_package"
    },
    error = function(e) {
      paste("render_failed:", conditionMessage(e))
    }
  )
} else {
  quarto_bin <- Sys.which("quarto")

  if (quarto_bin != "") {
    render_status <- tryCatch(
      {
        system2(
          quarto_bin,
          args = c("render", report_qmd, "--output", basename(report_html)),
          stdout = TRUE,
          stderr = TRUE
        )
        "rendered_with_quarto_cli"
      },
      error = function(e) {
        paste("render_failed:", conditionMessage(e))
      }
    )
  } else {
    render_status <- "quarto_not_available"
  }
}

cat("Reproducible proteomics report build complete.\n")
cat("\n")
cat("Report source:", report_qmd, "\n")
cat("Expected HTML:", report_html, "\n")
cat("Render status:", render_status, "\n")

if (render_status == "quarto_not_available") {
  cat("\n")
  cat("Quarto was not detected. Render manually with:\n")
  cat("quarto render", report_qmd, "\n")
}
