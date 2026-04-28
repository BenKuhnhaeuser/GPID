#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)

if (!length(args) %in% c(2, 3)) {
  stop("Usage: validation_preparations.R <validation_blast.tsv> <prepared_output.rds> [species_groups.csv]", call. = FALSE)
}

blast_file <- args[[1]]
prepared_output <- args[[2]]
species_groups_file <- if (length(args) == 3) args[[3]] else ""

required_columns <- c("gene", "query", "target", "pident", "length", "mismatch", "gapopen", "evalue", "bitscore")

suppressPackageStartupMessages({
  library(dplyr)
})

ids <- read.csv(blast_file, sep = "\t", stringsAsFactors = FALSE, check.names = FALSE)

missing_columns <- setdiff(required_columns, names(ids))
if (length(missing_columns) > 0) {
  stop(paste("BLAST file is missing required columns:", paste(missing_columns, collapse = ", ")), call. = FALSE)
}

if (nrow(ids) == 0) {
  stop("BLAST file does not contain any validation matches.", call. = FALSE)
}

species_name <- function(x) {
  parts <- strsplit(as.character(x), "_", fixed = TRUE)
  vapply(parts, function(part) {
    if (length(part) < 2) {
      NA_character_
    } else {
      paste(part[1:2], collapse = "_")
    }
  }, character(1))
}

genus_name <- function(x) {
  parts <- strsplit(as.character(x), "_", fixed = TRUE)
  vapply(parts, function(part) {
    if (length(part) < 1) {
      NA_character_
    } else {
      part[[1]]
    }
  }, character(1))
}

ids <- ids %>%
  mutate(
    query_sp = species_name(query),
    target_sp = species_name(target),
    query_samples = n_distinct(query)
  )

if (nzchar(species_groups_file)) {
  species_groups <- read.csv(species_groups_file, stringsAsFactors = FALSE, colClasses = "character")

  if (!identical(names(species_groups), c("genus_species", "species_group"))) {
    stop("Species groups file must use the header 'genus_species,species_group'.", call. = FALSE)
  }

  if (anyDuplicated(species_groups$genus_species) > 0) {
    duplicated_species <- unique(species_groups$genus_species[duplicated(species_groups$genus_species)])
    stop(
      "Species groups file contains duplicated species names: ",
      paste(duplicated_species, collapse = ", "),
      call. = FALSE
    )
  }

  ids <- left_join(ids, species_groups, by = c("query_sp" = "genus_species")) %>%
    rename(query_group = species_group)
  ids <- left_join(ids, species_groups, by = c("target_sp" = "genus_species")) %>%
    rename(target_group = species_group)
} else {
  ids <- ids %>%
    mutate(
      query_group = genus_name(query_sp),
      target_group = genus_name(target_sp)
    )
}

ids <- ids %>%
  mutate(
    id_correct = query_sp == target_sp,
    id_correct_group = query_group == target_group,
    id_correct_close = ifelse(
      id_correct,
      "correct",
      ifelse(is.na(id_correct_group), "wrong", ifelse(id_correct_group, "close", "wrong"))
    )
  )

factor_columns <- c(
  "gene", "query", "target", "query_sp", "target_sp",
  "query_group", "target_group", "id_correct", "id_correct_group",
  "id_correct_close"
)
ids[factor_columns] <- lapply(ids[factor_columns], factor)
ids$id_correct_close <- factor(ids$id_correct_close, c("correct", "close", "wrong"))

dir.create(dirname(prepared_output), recursive = TRUE, showWarnings = FALSE)
saveRDS(ids, prepared_output)

message("Prepared ", nrow(ids), " validation BLAST match(es) across ", length(unique(ids$query)), " sample(s).")
message("Output file written:")
message(prepared_output)
