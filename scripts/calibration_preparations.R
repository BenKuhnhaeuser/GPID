#!/usr/bin/env Rscript
# --- Purpose ---
# Prepare calibration RDS data by deriving species identities and match correctness from BLAST.

# --- Read CLI paths and expected BLAST columns ---
args <- commandArgs(trailingOnly = TRUE)

if (length(args) != 2) {
  stop("Usage: calibration_preparations.R <calibration_blast.tsv> <prepared_output.rds>", call. = FALSE)
}

blast_file <- args[[1]]
prepared_output <- args[[2]]

required_columns <- c("gene", "query", "target", "pident", "length", "mismatch", "gapopen", "evalue", "bitscore")

# --- Load BLAST matches ---
ids <- read.csv(blast_file, sep = "\t", stringsAsFactors = FALSE, check.names = FALSE)

missing_columns <- setdiff(required_columns, names(ids))
if (length(missing_columns) > 0) {
  stop(paste("BLAST file is missing required columns:", paste(missing_columns, collapse = ", ")), call. = FALSE)
}

# species_name(): Extract the first two underscore-delimited header fields as the species name.
species_name <- function(x) {
  parts <- strsplit(as.character(x), "_", fixed = TRUE)
  # Callback: Join the genus and species tokens, or return NA for incomplete headers.
  vapply(parts, function(part) {
    if (length(part) < 2) {
      NA_character_
    } else {
      paste(part[1:2], collapse = "_")
    }
  }, character(1))
}

# --- Derive exact-species correctness for calibration ---
# Calibration uses correct/wrong labels; group-based close labels are assigned
# later in the separate confidence workflow, not in threshold calibration.
ids$query_sp <- species_name(ids$query)
ids$target_sp <- species_name(ids$target)
ids$query_samples <- length(unique(ids$query))
ids$query_group <- NA_character_
ids$target_group <- NA_character_
ids$id_correct <- ids$query_sp == ids$target_sp
ids$id_correct_group <- NA
ids$id_correct_close <- ifelse(ids$id_correct, "correct", "wrong")

factor_columns <- c(
  "gene", "query", "target", "query_sp", "target_sp",
  "query_group", "target_group", "id_correct", "id_correct_group",
  "id_correct_close"
)
ids[factor_columns] <- lapply(ids[factor_columns], factor)
ids$id_correct_close <- factor(ids$id_correct_close, c("correct", "close", "wrong"))

# --- Save the prepared dataset ---
dir.create(dirname(prepared_output), recursive = TRUE, showWarnings = FALSE)
saveRDS(ids, prepared_output)

message("Prepared ", nrow(ids), " calibration BLAST match(es) across ", length(unique(ids$query)), " sample(s).")
message("Output file written:")
message(prepared_output)
