#!/usr/bin/env Rscript

# 08-prepare-functional-enrichment-inputs.R
#
# Prepare gene and protein identifier lists for functional enrichment.
#
# Usage:
#   Rscript scripts/R/08-prepare-functional-enrichment-inputs.R results/cleaned-protein-identifiers.tsv results

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
      "Rscript scripts/R/08-prepare-functional-enrichment-inputs.R <cleaned_protein_identifiers.tsv> [output_dir]",
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

all_genes_file <- file.path(output_dir, "enrichment-input-all-genes.tsv")
up_genes_file <- file.path(output_dir, "enrichment-input-upregulated-genes.tsv")
down_genes_file <- file.path(output_dir, "enrichment-input-downregulated-genes.tsv")

all_proteins_file <- file.path(output_dir, "enrichment-input-all-proteins.tsv")
up_proteins_file <- file.path(output_dir, "enrichment-input-upregulated-proteins.tsv")
down_proteins_file <- file.path(output_dir, "enrichment-input-downregulated-proteins.tsv")

summary_file <- file.path(output_dir, "enrichment-input-summary.tsv")

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

make_identifier_list <- function(tbl, id_col, output_col) {
  if (is.na(id_col)) {
    return(tibble(!!output_col := character()))
  }

  tbl %>%
    transmute(identifier = stringr::str_trim(as.character(.data[[id_col]]))) %>%
    filter(!is.na(identifier), identifier != "") %>%
    distinct(identifier) %>%
    arrange(identifier) %>%
    rename(!!output_col := identifier)
}

delimiter <- detect_delimiter(input_file)

cleaned_tbl <- readr::read_delim(
  file = input_file,
  delim = delimiter,
  show_col_types = FALSE,
  progress = FALSE,
  name_repair = "minimal"
)

column_names <- names(cleaned_tbl)

gene_col <- find_column(
  column_names,
  c("gene_symbol_clean", "gene_symbol", "gene_name", "Gene_Name", "symbol"),
  "gene symbol",
  required = FALSE
)

protein_col <- find_column(
  column_names,
  c("primary_protein_id", "protein_id_clean", "protein_id", "accession", "uniprot_id"),
  "protein identifier",
  required = FALSE
)

regulation_col <- find_column(
  column_names,
  c("regulation", "direction", "change_direction"),
  "regulation direction",
  required = FALSE
)

status_col <- find_column(
  column_names,
  c("identifier_status", "status"),
  "identifier status",
  required = FALSE
)

usable_tbl <- cleaned_tbl

if (!is.na(status_col)) {
  usable_tbl <- usable_tbl %>%
    filter(.data[[status_col]] %in% c("clean_identifier", "multiple_identifiers"))
}

if (!is.na(regulation_col)) {
  up_tbl <- usable_tbl %>%
    filter(.data[[regulation_col]] == "upregulated")

  down_tbl <- usable_tbl %>%
    filter(.data[[regulation_col]] == "downregulated")
} else {
  up_tbl <- usable_tbl[0, , drop = FALSE]
  down_tbl <- usable_tbl[0, , drop = FALSE]
}

all_genes <- make_identifier_list(usable_tbl, gene_col, "gene_symbol")
up_genes <- make_identifier_list(up_tbl, gene_col, "gene_symbol")
down_genes <- make_identifier_list(down_tbl, gene_col, "gene_symbol")

all_proteins <- make_identifier_list(usable_tbl, protein_col, "protein_id")
up_proteins <- make_identifier_list(up_tbl, protein_col, "protein_id")
down_proteins <- make_identifier_list(down_tbl, protein_col, "protein_id")

summary_tbl <- tibble(
  metric = c(
    "input_file",
    "total_cleaned_rows",
    "usable_rows",
    "gene_symbol_column",
    "protein_identifier_column",
    "regulation_column",
    "all_gene_identifiers",
    "upregulated_gene_identifiers",
    "downregulated_gene_identifiers",
    "all_protein_identifiers",
    "upregulated_protein_identifiers",
    "downregulated_protein_identifiers"
  ),
  value = c(
    input_file,
    as.character(nrow(cleaned_tbl)),
    as.character(nrow(usable_tbl)),
    ifelse(is.na(gene_col), "not_detected", gene_col),
    ifelse(is.na(protein_col), "not_detected", protein_col),
    ifelse(is.na(regulation_col), "not_detected", regulation_col),
    as.character(nrow(all_genes)),
    as.character(nrow(up_genes)),
    as.character(nrow(down_genes)),
    as.character(nrow(all_proteins)),
    as.character(nrow(up_proteins)),
    as.character(nrow(down_proteins))
  )
)

readr::write_tsv(all_genes, all_genes_file)
readr::write_tsv(up_genes, up_genes_file)
readr::write_tsv(down_genes, down_genes_file)

readr::write_tsv(all_proteins, all_proteins_file)
readr::write_tsv(up_proteins, up_proteins_file)
readr::write_tsv(down_proteins, down_proteins_file)

readr::write_tsv(summary_tbl, summary_file)

cat("Functional enrichment input preparation complete.\n")
cat("\n")
cat("Input file:", input_file, "\n")
cat("Usable rows:", nrow(usable_tbl), "\n")
cat("Gene symbol column:", ifelse(is.na(gene_col), "not detected", gene_col), "\n")
cat("Protein identifier column:", ifelse(is.na(protein_col), "not detected", protein_col), "\n")
cat("\n")
cat("Output files:\n")
cat(" -", all_genes_file, "\n")
cat(" -", up_genes_file, "\n")
cat(" -", down_genes_file, "\n")
cat(" -", all_proteins_file, "\n")
cat(" -", up_proteins_file, "\n")
cat(" -", down_proteins_file, "\n")
cat(" -", summary_file, "\n")
