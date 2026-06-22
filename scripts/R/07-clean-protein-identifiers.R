#!/usr/bin/env Rscript

# 07-clean-protein-identifiers.R
#
# Clean protein identifiers from ranked significant proteomics results.
#
# Usage:
#   Rscript scripts/R/07-clean-protein-identifiers.R results/ranked-significant-proteins.tsv results

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
      "Rscript scripts/R/07-clean-protein-identifiers.R <ranked_significant_proteins.tsv> [output_dir]",
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

cleaned_file <- file.path(output_dir, "cleaned-protein-identifiers.tsv")
protein_list_file <- file.path(output_dir, "annotation-ready-protein-list.tsv")
gene_list_file <- file.path(output_dir, "annotation-ready-gene-list.tsv")
string_input_file <- file.path(output_dir, "string-network-input.tsv")
summary_file <- file.path(output_dir, "identifier-cleaning-summary.tsv")

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

clean_accession <- function(x) {
  x <- as.character(x)
  x <- stringr::str_trim(x)

  # Remove common UniProt pipe format: sp|P12345|NAME or tr|Q12345|NAME
  x <- ifelse(
    stringr::str_detect(x, "^(sp|tr)\\|[^|]+\\|"),
    stringr::str_replace(x, "^(sp|tr)\\|([^|]+)\\|.*$", "\\2"),
    x
  )

  x <- stringr::str_replace_all(x, "^CON__", "")
  x <- stringr::str_replace_all(x, "^REV__", "")
  x <- stringr::str_replace_all(x, "^Reverse_", "")
  x <- stringr::str_trim(x)

  x
}

get_primary_identifier <- function(x) {
  x <- as.character(x)
  x <- stringr::str_split(x, pattern = "[;,| ]+", simplify = TRUE)[, 1]
  stringr::str_trim(x)
}

delimiter <- detect_delimiter(input_file)

ranked_tbl <- readr::read_delim(
  file = input_file,
  delim = delimiter,
  show_col_types = FALSE,
  progress = FALSE,
  name_repair = "minimal"
)

column_names <- names(ranked_tbl)

protein_id_col <- find_column(
  column_names,
  c("protein_id", "protein_id_original", "Protein_ID", "accession", "Accession", "uniprot_id", "UniProt_ID"),
  "protein identifier"
)

gene_symbol_col <- find_column(
  column_names,
  c("gene_symbol", "gene_name", "Gene_Name", "genes", "Gene", "symbol"),
  "gene symbol",
  required = FALSE
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

cleaned_tbl <- ranked_tbl %>%
  mutate(
    protein_id_original = as.character(.data[[protein_id_col]]),
    protein_id_clean = clean_accession(protein_id_original),
    primary_protein_id = get_primary_identifier(protein_id_clean),
    is_multi_identifier = stringr::str_detect(protein_id_original, "[;,]"),
    is_contaminant = stringr::str_detect(
      stringr::str_to_lower(protein_id_original),
      "con__|contaminant|keratin|trypsin"
    ),
    is_reverse = stringr::str_detect(
      stringr::str_to_lower(protein_id_original),
      "rev__|reverse|decoy"
    ),
    identifier_status = case_when(
      is.na(primary_protein_id) | primary_protein_id == "" ~ "missing_identifier",
      is_contaminant ~ "possible_contaminant",
      is_reverse ~ "possible_reverse_or_decoy",
      is_multi_identifier ~ "multiple_identifiers",
      TRUE ~ "clean_identifier"
    )
  )

if (!is.na(gene_symbol_col)) {
  cleaned_tbl <- cleaned_tbl %>%
    mutate(gene_symbol_clean = stringr::str_trim(as.character(.data[[gene_symbol_col]])))
} else {
  cleaned_tbl <- cleaned_tbl %>%
    mutate(gene_symbol_clean = NA_character_)
}

annotation_candidates <- cleaned_tbl %>%
  filter(
    identifier_status %in% c("clean_identifier", "multiple_identifiers"),
    !is.na(primary_protein_id),
    primary_protein_id != ""
  )

protein_list <- annotation_candidates %>%
  distinct(primary_protein_id) %>%
  arrange(primary_protein_id) %>%
  rename(protein_id = primary_protein_id)

gene_list <- annotation_candidates %>%
  filter(!is.na(gene_symbol_clean), gene_symbol_clean != "") %>%
  distinct(gene_symbol_clean) %>%
  arrange(gene_symbol_clean) %>%
  rename(gene_symbol = gene_symbol_clean)

string_input <- annotation_candidates %>%
  mutate(
    string_identifier = case_when(
      !is.na(gene_symbol_clean) & gene_symbol_clean != "" ~ gene_symbol_clean,
      !is.na(primary_protein_id) & primary_protein_id != "" ~ primary_protein_id,
      TRUE ~ protein_id_clean
    )
  )

if (!is.na(rank_col)) {
  string_input <- string_input %>% arrange(.data[[rank_col]])
}

string_keep_cols <- c(
  "string_identifier",
  "primary_protein_id",
  "gene_symbol_clean",
  "protein_id_original",
  "protein_id_clean",
  "identifier_status"
)

if (!is.na(rank_col)) {
  string_keep_cols <- c(rank_col, string_keep_cols)
}

if (!is.na(regulation_col)) {
  string_keep_cols <- c(string_keep_cols, regulation_col)
}

string_input <- string_input %>%
  select(any_of(string_keep_cols)) %>%
  distinct()

summary_tbl <- tibble(
  metric = c(
    "input_file",
    "input_rows",
    "unique_cleaned_protein_identifiers",
    "gene_symbols_available",
    "multi_identifier_rows",
    "possible_contaminant_rows",
    "possible_reverse_or_decoy_rows",
    "annotation_ready_protein_ids",
    "annotation_ready_gene_symbols"
  ),
  value = c(
    input_file,
    as.character(nrow(ranked_tbl)),
    as.character(n_distinct(cleaned_tbl$primary_protein_id, na.rm = TRUE)),
    as.character(sum(!is.na(cleaned_tbl$gene_symbol_clean) & cleaned_tbl$gene_symbol_clean != "")),
    as.character(sum(cleaned_tbl$is_multi_identifier, na.rm = TRUE)),
    as.character(sum(cleaned_tbl$is_contaminant, na.rm = TRUE)),
    as.character(sum(cleaned_tbl$is_reverse, na.rm = TRUE)),
    as.character(nrow(protein_list)),
    as.character(nrow(gene_list))
  )
)

readr::write_tsv(cleaned_tbl, cleaned_file)
readr::write_tsv(protein_list, protein_list_file)
readr::write_tsv(gene_list, gene_list_file)
readr::write_tsv(string_input, string_input_file)
readr::write_tsv(summary_tbl, summary_file)

cat("Protein identifier cleaning complete.\n")
cat("\n")
cat("Input file:", input_file, "\n")
cat("Input rows:", nrow(ranked_tbl), "\n")
cat("Protein ID column:", protein_id_col, "\n")
cat("Gene symbol column:", ifelse(is.na(gene_symbol_col), "not detected", gene_symbol_col), "\n")
cat("\n")
cat("Output files:\n")
cat(" -", cleaned_file, "\n")
cat(" -", protein_list_file, "\n")
cat(" -", gene_list_file, "\n")
cat(" -", string_input_file, "\n")
cat(" -", summary_file, "\n")
